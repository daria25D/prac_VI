#include <iostream>
#include <complex>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <string.h>
#include <omp.h>
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

complexd * generate(int n) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    long long size_array = power(2, n) / size;
    complexd * a = new complexd[size_array];
    double length = 0;
    long long i;
    unsigned int seed = time(0);
    //cout << seed << endl;
    seed *= rank + 1;
    //do omp for
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

double measure(complexd * noise, complexd * ideal, int n) {
    complexd M(0.0, 0.0);
    for (int i = 0; i < n; i++)
        M += conj(noise[i]) * ideal[i];
    return abs(M);
}

double normal_dis_gen() { 
    double S = 0.;
    for (int i = 0; i < 12; i++)  { 
        S += (double)rand()/RAND_MAX; 
    }
    return S-6.;
}

int main(int argc, char ** argv) {
    // ./task3 n eps iter_count threads output_file F_file <input_file>
    // if iter_count == 1 - only for time computation
    int n, k, iter_count, threads;
    double eps;
    bool out_flag = false, gen_flag = true;
    if (argc != 7 && argc != 8) {
        cout << "Not enough arguments\n";
        return -1;
    }
    if (strcmp(argv[5], "no") != 0) { //output mode
        out_flag = true;
    } 
    if (argc == 8) gen_flag = false;
    n = atoi(argv[1]);  
    eps = atof(argv[2]);
    iter_count = atoi(argv[3]);
    threads = atoi(argv[4]);
	omp_set_num_threads(threads);
    complexd ** U = new complexd*[2];
    U[0] = new complexd[2];
    U[1] = new complexd[2];
    U[0][0] = 1 / sqrt(2);
    U[0][1] =  U[0][0];
    U[1][0] =  U[0][0];
    U[1][1] = -U[0][0];
    srand(time(NULL));
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    long long size_array;
    int complex_size, int_size;
    MPI_Type_size(MPI_CXX_DOUBLE_COMPLEX, &complex_size);
    MPI_Type_size(MPI_INT, &int_size);
    double timer = 0.0;
    complexd * a, * b_noise, * b_ideal, * a_file, * a_ideal;
    size_array = power(2, n) / size;
    if (!gen_flag) {
        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, argv[7], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        MPI_File_read_all(fh, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_set_view(fh, rank * size_array * complex_size + int_size, MPI_CXX_DOUBLE_COMPLEX, 
                          MPI_CXX_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
        a_file = new complexd[size_array];
        a_file = normalize(a_file, size_array);
        MPI_File_read(fh, a_file, size_array, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        normalize(a_file, size_array);
        MPI_File_close(&fh);
    }
    double F;
    ofstream f;
    complexd ** H_e = new complexd*[2];
    H_e[0] = new complexd[2];
    H_e[1] = new complexd[2];
    if (iter_count > 1 && rank == 0) f.open(argv[6], ofstream::out);
    for (int i = 0; i < iter_count; i++) {
        if (gen_flag) a = generate(n);
        else {
            a = new complexd[size_array];
            memmove(a, a_file, sizeof(complexd) * size_array);
        }
        k = 1;
        a_ideal = new complexd[size_array];  
        memmove(a_ideal, a, sizeof(complexd) * size_array);

        for (int j = 0; j < n; j++) {
            if (rank == 0) { //generate H_e
                double xi = normal_dis_gen();
                double theta = xi * eps;
                H_e[0][0] =  (cos(theta) - sin(theta))/sqrt(2);
                H_e[0][1] =  (sin(theta) + cos(theta))/sqrt(2);
                H_e[1][0] =  H_e[0][1];
                H_e[1][1] = -H_e[0][0];
            }
            MPI_Bcast(H_e[0], 2, MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(H_e[1], 2, MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            double timer2 = MPI_Wtime();
            b_noise = quantum(a, n, H_e, k);
            timer2 = MPI_Wtime() - timer2;
            timer += timer2;
            memmove(a, b_noise, sizeof(complexd) * size_array);
            if (iter_count > 1) {
                b_ideal = quantum(a_ideal, n, U, k);
                memmove(a_ideal, b_ideal, sizeof(complexd) * size_array);
            }
            if (k != n) {
                delete[] b_ideal;
                delete[] b_noise;
            } 
            k++;               
        }
        if (iter_count > 1) {
            F = 0.0;
            double F_local = measure(b_noise, b_ideal, size_array);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Reduce(&F_local, &F, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if (rank == 0) {
                f << 1 - F << endl;
            }  
            delete[] b_ideal;
        }
        delete[] a;
        delete[] a_ideal;
        delete[] b_noise;
    }
   
    if (out_flag && iter_count == 1) { //output mode, only if iter_count == 1
        //argv[5] - output file
        MPI_File out;
        MPI_File_open(MPI_COMM_WORLD, argv[5], MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &out);

        if (rank == 0) {
            MPI_File_write(out, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
            MPI_File_write(out, &k, 1, MPI_INT, MPI_STATUS_IGNORE);
        }
        MPI_File_set_view(out, rank * size_array * complex_size + 2 * int_size, MPI_CXX_DOUBLE_COMPLEX, 
                          MPI_CXX_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
        MPI_File_write(out, b_noise, size_array, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        MPI_File_close(&out);
    } 
    double timer2;
    MPI_Reduce(&timer, &timer2, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (rank == 0 && iter_count == 1) cout << timer2 << endl;

    if (iter_count > 1) {
        if (rank == 0) f.close();
    }
    if (!gen_flag) delete[] a_file;
    MPI_Finalize();
    return 0;
}