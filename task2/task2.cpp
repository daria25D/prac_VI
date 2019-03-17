#include <iostream>
#include <complex>
#include <cstdlib>
#include <random>
#include <cmath>	
#include <ctime>
#include <fstream>
#include <string.h>
#include "mpi.h"

typedef std::complex<double> complexd;
using namespace std;

long long power(int base, int deg) {
	long long res = 1;
	if (base == 2) {
		for (int i = 0; i < deg; i++)
			res <<= 1;
	} else {
		for (int i = 0; i < deg; i++)
			res *= base;
	}
	return res;
}
complexd * quantum(complexd * a, int n, complexd ** U, int k) {

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
    bool flag_exchage = false;
    if (pow(2, k) <= size) {
        flag_exchage = true;
    }
    if (flag_exchage) {          		
        idx1 = (rank * size_of_a) & n0;
		idx2 = (rank * size_of_a) | n1;
        for (long long j = 0; j < size; j++) {
            if (idx1 >= j * size_of_a && idx1 < (j + 1) * size_of_a) 
                rank1 = j;
            if (idx2 >= j * size_of_a && idx2 < (j + 1) * size_of_a) 
                rank2 = j;
            if (rank1 != 0 && rank2 != 0) break;
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
	for (long long i = 0; i < size_of_a; i++) {
		long long idx_ik = ((i + rank * size_of_a) >> (n - k)) % 2;
		idx1 = (i + rank * size_of_a) & n0;
		idx2 = (i + rank * size_of_a) | n1;

        if (flag_exchage) {
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
    for (i = 0; i < size_array; i++) {
        a[i].real(rand_r(&seed));
        a[i].imag(rand_r(&seed));
        length += abs(a[i]) * abs(a[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&length, &length, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    length = sqrt(length);
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


int main(int argc, char ** argv) {
    // ./task2 n k gen output
    // or
    // ./task2 file_in file_out
	int n, k;
    bool gen_flag = false;
    if (argc != 3 && argc != 5) {
        cout << "Not enough arguments\n";
        return -1;
    }
    if (argc == 5 && strcmp(argv[3], "gen") == 0) { //generation mod
        gen_flag = true;
	    n = atoi(argv[1]);
	    k = atoi(argv[2]);
    }   
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	long long size_array;
    int complex_size, int_size;
    MPI_Type_size(MPI_CXX_DOUBLE_COMPLEX, &complex_size);
    MPI_Type_size(MPI_INT, &int_size);
    double timer1 = 0.0;
    complexd * a;
	if (gen_flag) { //generation mode
        size_array = power(2, n) / size;
        timer1 = MPI_Wtime();
	    a = generate(n);
        timer1 = MPI_Wtime() - timer1;
    } else {
        //read from file
        MPI_File f;
        MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
        MPI_File_read_all(f, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_read_all(f, &k, 1, MPI_INT, MPI_STATUS_IGNORE);
        size_array = power(2, n) / size;
        MPI_File_set_view(f, rank * size_array * complex_size + 2 * int_size, MPI_CXX_DOUBLE_COMPLEX, 
                          MPI_CXX_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
        a = new complexd[size_array];
        MPI_File_read(f, a, size_array, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        normalize(a, size_array);
        MPI_File_close(&f);
    }
	complexd ** U = new complexd*[2];
	U[0] = new complexd[2];
	U[1] = new complexd[2];
	U[0][0] = 1 / sqrt(2);
	U[0][1] = 1 / sqrt(2);
	U[1][0] = 1 / sqrt(2);
	U[1][1] = -1 / sqrt(2);

	double timer2 = MPI_Wtime();
	complexd * b = quantum(a, n, U, k);
	timer2 = MPI_Wtime() - timer2;

    if (gen_flag && strcmp(argv[4], "no") != 0) { //generation mode, output mode
        //argv[4] - output file
        MPI_File out;
        MPI_File_open(MPI_COMM_WORLD, argv[4], MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &out);
        if (rank == 0) {
            MPI_File_write(out, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
            MPI_File_write(out, &k, 1, MPI_INT, MPI_STATUS_IGNORE);
        }
        MPI_File_set_view(out, rank * size_array * complex_size + 2 * int_size, MPI_CXX_DOUBLE_COMPLEX, 
                          MPI_CXX_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
        MPI_File_write(out, b, size_array, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        MPI_File_close(&out);

    } else if (!gen_flag && strcmp(argv[2], "no") != 0) { //non-generation mode, output mode
        //argv[2] - output file
        MPI_File out;
        MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &out);
        if (rank == 0) {
            MPI_File_write(out, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
            MPI_File_write(out, &k, 1, MPI_INT, MPI_STATUS_IGNORE);
        }
        MPI_File_set_view(out, rank * size_array * complex_size + 2 * int_size, MPI_CXX_DOUBLE_COMPLEX, 
                          MPI_CXX_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
        MPI_File_write(out, b, size_array, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        MPI_File_close(&out);
    } else {
        //no output
    }
	if (rank == 0) cout << timer2 << endl;
	delete a;
	delete b;
    MPI_Finalize();
	return 0;
}