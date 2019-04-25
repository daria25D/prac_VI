#pragma once
#ifndef TASK4_QUANTUM_H_
#define TASK4_QUANTUM_H_

#include <omp.h>
#include <iostream>
#include <complex>
#include <cstdlib>
#include <random>
#include <cmath>
#include <ctime>
#include "./mpi.h"

typedef std::complex<double> complexd;
using namespace std;

const complexd U_[2][2] = {{1 / sqrt(2), 1 / sqrt(2)}, {1 / sqrt(2), -1 / sqrt(2)}};

uint64_t power(uint64_t base, uint64_t deg) {
    uint64_t res = 1;
    if (base == 2) {
        for (uint64_t i = 0; i < deg; i++)
            res <<= 1;
    } else {
        for (uint64_t i = 0; i < deg; i++)
            res *= base;
    }
    return res;
}

complexd * normalize(complexd * a, uint64_t size_array) {
    double length = 0.0;
    for (uint64_t i = 0; i < size_array; i++) {
        length += abs(a[i]) * abs(a[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&length, &length, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    length = sqrt(length);
    for (uint64_t i = 0; i < size_array; i++) {
        a[i] /= length;
    }
	return a;
}

complexd * generate(int n) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    uint64_t size_array = power(2, n) / size;
    complexd * a = new complexd[size_array];
    double length = 0;
    uint64_t i;
    unsigned int seed = time(0);
    //cout << size << endl;
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
    for (i = 0; i < size_array; i++) {
        a[i] /= length;
    }
    return a;
}

complexd * quantum(complexd * a, int n, const complexd U[][2], int k) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    uint64_t size_of_a = power(2, n);
    if (size != 1) size_of_a /= size;
    uint64_t n1 = 1;
    uint64_t idx1, idx2;        
    int rank1 = 0, rank2 = 0;
    for (int j = 0; j < n - k; j++)
        n1 <<= 1;
    uint64_t n0 = 1;
    for (int j = 0; j < n; j++) {
        if (j == k - 1)
            n0 <<= 1;
        else
            n0 = (n0 << 1) | 1;
    }
    complexd * b = new complexd[size_of_a];
    complexd * a_swap = new complexd;
    bool flag_exchange = false;
    if (pow(2, k) <= size) {
        flag_exchange = true;
    }
    if (flag_exchange) {          
        idx1 = (rank * size_of_a) & n0;
        idx2 = (rank * size_of_a) | n1;
        rank1 = idx1/static_cast<int>(size_of_a);
        rank2 = idx2/static_cast<int>(size_of_a);
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
    for (uint64_t i = 0; i < size_of_a; i++) {
        uint64_t idx_ik = ((i + rank * size_of_a) >> (n - k)) % 2;
        idx1 = (i + rank * size_of_a) & n0;
        idx2 = (i + rank * size_of_a) | n1;
        if (flag_exchange) {
            if (rank == rank1) {
                b[i] = U[idx_ik][0] * a[idx1 - rank * size_of_a] + U[idx_ik][1] * a_swap[idx2 - rank2 * size_of_a];
            } else if (rank == rank2) {
                b[i] = U[idx_ik][0] * a_swap[idx1 - rank1 * size_of_a] + U[idx_ik][1] * a[idx2 - rank * size_of_a];    
            }                
        } else { // no messages needed
            b[i] = U[idx_ik][0] * a[idx1 - rank * size_of_a] + U[idx_ik][1] * a[idx2 - rank * size_of_a]; 
        }
    }
    return b;
}

complexd * quantum4x4(complexd * a, int n, complexd U[][4], int k, int l) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    uint64_t size_of_a = power(2, n);
    if (size != 1) size_of_a /= size;
    uint64_t n11 = 1, n10 = 1;
    uint64_t idx1, idx2, idx3, idx4;        
    int rank1 = 0, rank2 = 0, rank3 = 0, rank4 = 0;
    if (l < k) {
        complexd swap[4];
        memcpy(&swap, &U[1], 4 * sizeof(complexd));
        memcpy(&U[1], &U[2], 4 * sizeof(complexd));
        memcpy(&U[2], &swap, 4 * sizeof(complexd));
        int sw = k;
        k = l;
        l = sw;
    }
    for (int j = 0; j < n - k; j++)
        n11 <<= 1;
    n10 = n11;
    uint64_t nl = 1l<<(n-l);
    n11 |= nl;
    uint64_t n00 = 1, n01;
    for (int j = 0; j < n; j++) {
        if (j == k - 1 || j == l - 1)
            n00 <<= 1;
        else
            n00 = (n00 << 1) | 1;
    }
    n01 = nl;
    //cout << n00 << " " << n01 << " " << n10 << " " << n11 << endl;
    complexd * b = new complexd[size_of_a];
    complexd * a_swap_l = new complexd;
    complexd * a_swap_1 = new complexd;
    complexd * a_swap_2 = new complexd;
    complexd * a_swap_3 = new complexd;
    bool flag_exchange = false, flag_l_only = false;
    if (pow(2, k) <= size || pow(2, l) <= size) {
        flag_exchange = true;
        // no messages - false
        // only for index l - true
        // for both k and l - true
    }
    if (flag_exchange && pow(2, k) > size) {
        flag_l_only = true;
    }
    if (flag_exchange) {          
        idx1 = (rank * size_of_a) & n00;
        idx2 = ((rank * size_of_a) & n00) | n01;
        idx3 = ((rank * size_of_a) & n00) | n10;
        idx4 = (rank * size_of_a) | n11;
        rank1 = idx1/static_cast<int>(size_of_a); // 00
        rank2 = idx2/static_cast<int>(size_of_a); // 01 different from rank1 <=> flag_l_only
        rank3 = idx3/static_cast<int>(size_of_a); // 10 different from rank1 <=> flag_exchange
        rank4 = idx4/static_cast<int>(size_of_a); // 11 different from rank3 <=> flag_l_only
                                                  // 11 different from rank2 <=> flag_exchange
                                                  // !flag_l_only <=> rank1 = rank2, rank3 = rank4
                                                  // !flag_exchange <=> rank1 = rank2 = rank3 = rank4
                                                  // flag_l_only && flag_exchange <=> rank1 != rank2 != rank3 != rank4
        if (flag_l_only) {
            a_swap_l = new complexd[size_of_a];
            if (rank == rank1) {
                MPI_Send(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank3, 0, MPI_COMM_WORLD);
                MPI_Recv(a_swap_l, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else if (rank == rank3) {
                MPI_Recv(a_swap_l, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank1, 0, MPI_COMM_WORLD);
            }
        } else if (flag_exchange && !flag_l_only) {
            a_swap_1 = new complexd[size_of_a];
            a_swap_2 = new complexd[size_of_a];
            a_swap_3 = new complexd[size_of_a];
            MPI_Request req;
            if (rank == rank1) {
                MPI_Isend(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank2, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Isend(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank3, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Isend(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank4, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Irecv(a_swap_1, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank4, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Irecv(a_swap_2, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank2, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Irecv(a_swap_3, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank3, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);                
            } else if (rank == rank2) {
                MPI_Isend(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank1, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Isend(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank3, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Isend(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank4, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Irecv(a_swap_1, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank1, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Irecv(a_swap_2, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank4, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Irecv(a_swap_3, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank3, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req); 
             } else if (rank == rank3) {
                MPI_Isend(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank1, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Isend(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank2, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Isend(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank4, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Irecv(a_swap_1, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank1, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Irecv(a_swap_2, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank2, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Irecv(a_swap_3, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank4, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req); 
             } else if (rank == rank4) {
                MPI_Isend(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank1, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Isend(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank2, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Isend(a, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank3, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Irecv(a_swap_1, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank1, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Irecv(a_swap_2, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank2, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                MPI_Irecv(a_swap_3, size_of_a, MPI_CXX_DOUBLE_COMPLEX, rank3, 0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req); 
             }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
// #pragma omp parallel shared(size_of_a, b, a_swap_l, rank, rank1, rank2, n00, n11, flag_exchange)
// #pragma omp for
    for (uint64_t i = 0; i < size_of_a; i++) {
        uint64_t idx_ik = ((((i + rank * size_of_a) >> (n - k)) % 2) << 1) | 
                           (((i + rank * size_of_a) >> (n - l)) % 2);
        idx4 = (i + rank * size_of_a) & n00;
        idx3 = ((i + rank * size_of_a) & n00) | n01;
        idx1 = ((i + rank * size_of_a) & n00) | n10;
        idx2 = (i + rank * size_of_a) | n11;

        // cout << idx_ik << " " << idx1 << " " << idx2 << " " << idx3 << " " << idx4 << " " << rank << endl;
        if (flag_exchange) {
            if (flag_l_only) {
                // cout << "flag_l_only" << endl;
                if (rank == rank1) {
                    b[i] = U[idx_ik][0] * a[idx1 - rank * size_of_a] + 
                           U[idx_ik][1] * a[idx2 - rank * size_of_a] +
                           U[idx_ik][2] * a_swap_l[idx3 - rank3 * size_of_a] + 
                           U[idx_ik][3] * a_swap_l[idx4 - rank3 * size_of_a];
                } else if (rank == rank3) {
                    b[i] = U[idx_ik][0] * a_swap_l[idx1 - rank1 * size_of_a] + 
                           U[idx_ik][1] * a_swap_l[idx2 - rank1 * size_of_a] + 
                           U[idx_ik][2] * a[idx3 - rank * size_of_a] +
                           U[idx_ik][3] * a[idx4 - rank * size_of_a];    
                }    
            } else {
                // cout << "flag_exchange" << endl;
                if (rank == rank1) {
                    b[i] = U[idx_ik][0] * a[idx1 - rank * size_of_a] +
                           U[idx_ik][1] * a_swap_2[idx2 - rank2 * size_of_a] + 
                           U[idx_ik][2] * a_swap_3[idx3 - rank3 * size_of_a] + 
                           U[idx_ik][3] * a_swap_1[idx4 - rank4 * size_of_a];
                } else if (rank == rank2) {
                    b[i] = U[idx_ik][0] * a_swap_1[idx1 - rank1 * size_of_a] +
                           U[idx_ik][1] * a[idx2 - rank * size_of_a] + 
                           U[idx_ik][2] * a_swap_3[idx3 - rank3 * size_of_a] + 
                           U[idx_ik][3] * a_swap_2[idx4 - rank4 * size_of_a];
                } else if (rank == rank3) {
                    b[i] = U[idx_ik][0] * a_swap_1[idx1 - rank1 * size_of_a] +
                           U[idx_ik][1] * a_swap_2[idx2 - rank2 * size_of_a] + 
                           U[idx_ik][2] * a[idx3 - rank * size_of_a] + 
                           U[idx_ik][3] * a_swap_3[idx4 - rank4 * size_of_a];
                } else if (rank == rank4) {
                    b[i] = U[idx_ik][0] * a_swap_1[idx1 - rank1 * size_of_a] +
                           U[idx_ik][1] * a_swap_2[idx2 - rank2 * size_of_a] + 
                           U[idx_ik][2] * a_swap_3[idx3 - rank3 * size_of_a] + 
                           U[idx_ik][3] * a[idx4 - rank * size_of_a];
                }
            }      
        } else { // no messages needed
            // cout << "no messages" << endl;
            b[i] = U[idx_ik][0] * a[idx1 - rank * size_of_a] +
                   U[idx_ik][1] * a[idx2 - rank * size_of_a] + 
                   U[idx_ik][2] * a[idx3 - rank * size_of_a] +
                   U[idx_ik][3] * a[idx4 - rank * size_of_a]; 
        }
    }
    return b;
}

complexd * quantum_Adamar(complexd * a, int n, int k) {
    return quantum(a, n, U_, k);
}

complexd * quantum_nAdamar(complexd * a, int n) {
    complexd * b = new complexd;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    uint64_t size_array = power(2, n) / size;
    for (int k = 0; k < n; k++) {
        b = quantum(a, n, U_, k);
        memmove(a, b, sizeof(complexd) * size_array);
        if (k != n) delete[] b;
    }
    return b;
}

complexd * quantum_Rotate(complexd * a, int n, int k, double phi) {
    const complexd R[2][2] = {{1, 0}, 
                              {0, complexd(cos(phi), sin(phi))}};
    return quantum(a, n, R, k);
}

complexd * quantum_Not(complexd * a, int n, int k) {
    const complexd N[2][2] = {{0, 1},
                              {1, 0}};
    return quantum(a, n, N, k);
}

complexd * quantum_CNot(complexd * a, int n, int k, int l) {
    complexd C[4][4] = {{1, 0, 0, 0}, 
                        {0, 1, 0, 0},
                        {0, 0, 0, 1}, 
                        {0, 0, 1, 0}};
    return quantum4x4(a, n, C, k, l);
}

complexd * quantum_CRotate(complexd * a, int n, int k, int l, double phi) {
    complexd R[4][4] = {{1, 0, 0, 0}, 
                        {0, 1, 0, 0}, 
                        {0, 0, 1, 0}, 
                        {0, 0, 0, complexd(cos(phi), sin(phi))}};
    return quantum4x4(a, n, R, k, l);
}

void File_MPI_Write(complexd * b, uint64_t size_array, int k, int l, char * filename, int rank) {
    MPI_File out;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &out);

    int complex_size, long_size;
    MPI_Type_size(MPI_CXX_DOUBLE_COMPLEX, &complex_size);
    MPI_Type_size(MPI_LONG_LONG, &long_size);
    int int_size;
    MPI_Type_size(MPI_INT, &int_size);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) { // printed single for nAdamar
        MPI_File_write(out, &size_array, 1, MPI_LONG_LONG, MPI_STATUS_IGNORE);
    }
    if (k != 0) { // printed for H, Rot, Not, CRot, CNot - only (or first) index
        if (rank == 0) MPI_File_write(out, &k, 1, MPI_INT, MPI_STATUS_IGNORE);
        long_size += int_size;
    }
    if (l != 0) { // printed for CRot, CNot - second index
        if (rank == 0) MPI_File_write(out, &l, 1, MPI_INT, MPI_STATUS_IGNORE);
        long_size += int_size;
    }
    MPI_File_set_view(out, rank * (size_array/size) * complex_size + long_size, MPI_CXX_DOUBLE_COMPLEX, 
                        MPI_CXX_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
    MPI_File_write(out, b, size_array/size, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    MPI_File_close(&out);

}

#endif // TASK4_QUANTUM_H_
