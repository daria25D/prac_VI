#include <iostream>
#include <complex>
#include <cstdlib>
#include <random>
#include <cmath>	
#include <ctime>
#include <fstream>
#include "mpi.h"

typedef std::complex<double> complexd;
using namespace std;

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
complexd * quantum(complexd * a, long long n, complexd ** U, long long k) {
    /**
    TODO
    messages between processes (send parts of array a)
    */
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	long long size_of_a = power(2, n) / size;
    long long n1 = 1;
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
	for (long long i = 0; i < size_of_a; i++) {
		long long idx1, idx2;
		long long idx_ik = ((i + rank * size_of_a) >> (n - k)) % 2;
		idx1 = (i + rank * size_of_a) & n0;
		idx2 = (i + rank * size_of_a) | n1;
        //send data
        //if size = 1  - no messages needed
        //else - need to compute which data to send
        if (size == 1)
		    b[i] = U[idx_ik][0] * a[idx1] + U[idx_ik][1] * a[idx2]; 
        else {
            int rank1 = 0, rank2 = 0;
            for (long long j = 0; j < size; j++) {
                if (idx1 >= j * size_of_a && idx1 < (j + 1) * size_of_a) 
                    rank1 = j;
                if (idx2 >= j * size_of_a && idx2 < (j + 1) * size_of_a) 
                    rank2 = j;
                if (rank1 != 0 && rank2 != 0) break;
            }
            complexd a_idx;
            double a_idx1_real, a_idx1_imag, a_idx2_real, a_idx2_imag;
            if (rank1 != rank) { //need to send a[idx2] and receive a[idx1]
                a_idx2_real = a[idx2 - rank * size_of_a].real();
                a_idx2_imag = a[idx2 - rank * size_of_a].imag();
                MPI_Status status;
                MPI_Send(&a_idx2_real, 1, MPI_DOUBLE, rank1, 0, MPI_COMM_WORLD);
                MPI_Recv(&a_idx1_real, 1, MPI_DOUBLE, rank1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&a_idx2_imag, 1, MPI_DOUBLE, rank1, 0, MPI_COMM_WORLD);
                MPI_Recv(&a_idx1_imag, 1, MPI_DOUBLE, rank1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                a_idx.real(a_idx1_real);
                a_idx.imag(a_idx1_imag);
                b[i] = U[idx_ik][0] * a_idx + U[idx_ik][1] * a[idx2 - rank * size_of_a];
            } else if (rank2 != rank) { //need to send a[idx1] and receive a[idx2]
                a_idx1_real = a[idx1 - rank * size_of_a].real();
                a_idx1_imag = a[idx1 - rank * size_of_a].imag();
                MPI_Status status;
                MPI_Recv(&a_idx2_real, 1, MPI_DOUBLE, rank2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&a_idx1_real, 1, MPI_DOUBLE, rank2, 0, MPI_COMM_WORLD);
                MPI_Recv(&a_idx2_imag, 1, MPI_DOUBLE, rank2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&a_idx1_imag, 1, MPI_DOUBLE, rank2, 0, MPI_COMM_WORLD);
                a_idx.real(a_idx2_real);
                a_idx.imag(a_idx2_imag);
                b[i] = U[idx_ik][0] * a[idx1 - rank * size_of_a] + U[idx_ik][1] * a_idx;
            } else { // no messages needed
		        b[i] = U[idx_ik][0] * a[idx1 - rank * size_of_a] + U[idx_ik][1] * a[idx2 - rank * size_of_a]; 
            }
        }
	}
	return b;
}

complexd * generate(long long n) {
    /*
    TODO
    partition of array a between processes

    UPD
    done, maybe needs corrrection
    */
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
    MPI_Allreduce(&length, &length, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    length = sqrt(length);
    for (i = 0; i < size_array; i++)
        a[i] /= length;
	return a;
}

int main(int argc, char ** argv) {
	long long n, k;
    bool gen_flag = false;
    if (argc == 4 && argv[3] == "gen") gen_flag = true;
	n = atoi(argv[1]);
	k = atoi(argv[2]);
	complexd ** U = new complexd*[2];
	U[0] = new complexd[2];
	U[1] = new complexd[2];
	U[0][0] = 1 / sqrt(2);
	U[0][1] = 1 / sqrt(2);
	U[1][0] = 1 / sqrt(2);
	U[1][1] = -1 / sqrt(2);
    MPI_Init(&argc, &argv);
    double timer1 = 0.0;
	if (gen_flag) timer1 = MPI_Wtime();
	complexd * a = generate(n);
    if (gen_flag) timer1 = MPI_Wtime() - timer1;
	double timer2 = MPI_Wtime();
	complexd * b = quantum(a, n, U, k);
	timer2 = MPI_Wtime() - timer2;
	//long long size = power(2, n);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	ofstream f("vectors.txt");
    /**
    TODO
    MPI file usage
    */
	// for (long long i = 0; i < size; i++) {
	// 	f << a[i] << "\n";
	// }
	// f << endl;
	// for (long long i = 0; i < size; i++) {
	// 	f << b[i] << "\n";
	// }

	if (rank == 0) cout << timer1 + timer2 << endl;
	delete a;
	delete b;
    MPI_Finalize();
	return 0;
}