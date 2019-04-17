#pragma once
#ifndef QUANTUM_H
#define QUANTUM_H

#include <iostream>
#include <complex>
#include <cstdlib>
#include <random>
#include <cmath>	
#include <ctime>
#include "mpi.h"
#include <omp.h>

typedef std::complex<double> complexd;
using namespace std;

const complexd U[2][2] = {{1 / sqrt(2), 1 / sqrt(2)}, {1 / sqrt(2), -1 / sqrt(2)}};

long long power(long long base, long long deg) {
	long long res = 1;
	if (base == 2) {
		for (long long i = 0; i < deg; i++)
			res <<= 1;
	} else {
		for (long long i = 0; i < deg; i++)
			res *= base;
	}
	return res;
}

complexd * generate(int n) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    long long size_array = power(2, n) / size;
    complexd * a = new complexd[size_array];
    double length = 0;
    long long i;
    unsigned int seed = time(0);
    seed *= rank + 1;
#pragma omp parallel shared(size, a, length) firstprivate(seed)
    {
        seed += omp_get_thread_num();
	#pragma omp for reduction(+:length) private(i)
        for (i = 0; i < size_array; i++) {
            a[i].real(rand_r(&seed));
            a[i].imag(rand_r(&seed));
            length += abs(a[i]) * abs(a[i]);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&length, &length, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#pragma omp single
    length = sqrt(length);
#pragma omp for private(i)
    for (i = 0; i < size_array; i++)
    a[i] /= length;
    return a;
}

complexd * normalize(complexd * a, long long size_array) {
    double length = 0.0;
    for (long long i = 0; i < size_array; i++) {
        length += abs(a[i]) * abs(a[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&length, &length, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    length = sqrt(length);
    for (long long i = 0; i < size_array; i++)
        a[i] /= length;
	return a;
}

complexd * quantum(complexd * a, int n, const complexd U[][2], int k) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    long long size_of_a = power(2, n);
    if (size != 1) size_of_a /= size;
    long long n1 = 1;
    long long idx1, idx2;        
    int rank1 = 0, rank2 = 0;
    for (long long j = 0; j < n - k; j++) 
        n1 <<= 1;
    long long n0 = 1;
    for (long long j = 0; j < n; j++) {
        if (j == k - 1)
            n0 <<= 1;
        else
            n0 = (n0 << 1) | 1;
    }
    complexd * b = new complexd[size_of_a];
    complexd * a_swap;
    bool flag_exchange = false;
    if (pow(2, k) <= size) {
        flag_exchange = true;
    }
    if (flag_exchange) {          
        idx1 = (rank * size_of_a) & n0;
        idx2 = (rank * size_of_a) | n1;
        bool flag = true;
#pragma omp parallel for shared(size_of_a, rank, rank1, rank2, size, idx1, idx2, flag)
        for (long long j = 0; j < size; j++) {
            if (flag) { 
                if (idx1 >= j * size_of_a && idx1 < (j + 1) * size_of_a) 
                    rank1 = j;
                if (idx2 >= j * size_of_a && idx2 < (j + 1) * size_of_a) 
                    rank2 = j;
                if (rank1 != 0 && rank2 != 0) flag = false;
            }
        }
        a_swap = new complexd[size_of_a];
        if (rank == rank1) {
            MPI_Send(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank2, 0, MPI_COMM_WORLD);
            MPI_Recv(a_swap, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else if (rank == rank2) {
            MPI_Recv(a_swap, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank1, 0, MPI_COMM_WORLD);
        }
    }
#pragma omp parallel shared(size_of_a, b, a_swap, rank, rank1, rank2, n0, n1, flag_exchange)
#pragma omp for
    for (long long i = 0; i < size_of_a; i++) {
        long long idx_ik = ((i + rank * size_of_a) >> (n - k)) % 2;
        idx1 = (i + rank * size_of_a) & n0;
        idx2 = (i + rank * size_of_a) | n1;
        if (flag_exchange) {
            if (rank == rank1) {
                b[i] = U[idx_ik][0] * a[idx1 - rank * size_of_a] + U[idx_ik][1] * a_swap[idx2 - rank2 * size_of_a];
            }
            else if (rank == rank2) {
                b[i] = U[idx_ik][0] * a_swap[idx1 - rank1 * size_of_a] + U[idx_ik][1] * a[idx2 - rank * size_of_a];    
            }                
        } else { // no messages needed
            b[i] = U[idx_ik][0] * a[idx1 - rank * size_of_a] + U[idx_ik][1] * a[idx2 - rank * size_of_a]; 
        }
    }
    return b;
}

complexd * quantum_Adamar(complexd * a, int n, int k) {
    return quantum(a, n, U, k);
}

complexd * quantum_nAdamar(complexd * a, int n) {
    complexd * b;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    long long size_array = power(2, n) / size;
    for (int k = 0; k < n; k++) {
        b = quantum(a, n, U, k);
        memmove(a, b, sizeof(complexd) * size_array);
        if (k != n) delete[] b;
    }
    return b;
}

complexd * quantum_Rotate(complexd * a, int n, int k, double phi) {
    const complexd R[2][2] = {{1, 0}, {0, complexd(cos(phi), sin(phi))}};
    return quantum(a, n, R, k);
}

complexd * quantum_Not(complexd * a, int n, int k) {
    const complexd N[2][2] = {{0, 1}, {1, 0}};
    return quantum(a, n, N, k);
}

void File_MPI_Write(complexd * b, long long size_array, int k, char * filename, int rank) {
    MPI_File out;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &out);

    int complex_size, long_size;
    MPI_Type_size(MPI_CXX_DOUBLE_COMPLEX, &complex_size);
    MPI_Type_size(MPI_LONG_LONG, &long_size);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        MPI_File_write(out, &size_array, 1, MPI_LONG_LONG, MPI_STATUS_IGNORE);
        if (k != 0) {
            MPI_File_write(out, &k, 1, MPI_LONG_LONG, MPI_STATUS_IGNORE);
            int int_size;
            MPI_Type_size(MPI_INT, &int_size);
            long_size += int_size;
        }
    }
    MPI_File_set_view(out, rank * (size_array/size) * complex_size + long_size, MPI_CXX_DOUBLE_COMPLEX, 
                        MPI_CXX_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
    MPI_File_write(out, b, size_array/size, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    MPI_File_close(&out);
}

#endif