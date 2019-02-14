#include <iostream>
#include <complex>
#include <cstdlib>
#include <random>
#include <cmath>	
#include <ctime>
#include <fstream>
#include <omp.h>

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
	long long size_of_a = power(2, n);
	complexd * b = new complexd[size_of_a];
	#pragma omp parallel shared(size_of_a, b)
	#pragma omp for 
	for (long long i = 0; i < size_of_a; i++) {
		long long idx1, idx2;
		long long idx_ik = (i >> (n - k)) % 2;
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
		idx1 = i & n0;
		idx2 = i | n1;
		b[i] = U[idx_ik][0] * a[idx1] + U[idx_ik][1] * a[idx2];
	}
	return b;
}

complexd * generate(long long n) {
	long long size = power(2, n);
	complexd * a = new complexd[size];
	double length = 0;
	long long i;
	unsigned int seed = time(0);
#pragma omp parallel shared(size, a, length)
#pragma omp for reduction(+:length) private(i, seed)
	for (i = 0; i < size; i++) {
		seed += omp_get_thread_num() + i + 1 + size/omp_get_num_threads();
		a[i].real(rand_r(&seed));
		a[i].imag(rand_r(&seed));
		length += abs(a[i]) * abs(a[i]);
	}
	length = sqrt(length);
#pragma omp for private(i)
	for (i = 0; i < size; i++)
		a[i] /= length;
	return a;
}

int main(int argc, char ** argv) {
	long long n, k;
	n = atoi(argv[1]);
	k = atoi(argv[2]);
	int n_threads = atoi(argv[3]);
	complexd ** U = new complexd*[2];
	U[0] = new complexd[2];
	U[1] = new complexd[2];
	U[0][0] = 1 / sqrt(2);
	U[0][1] = 1 / sqrt(2);
	U[1][0] = 1 / sqrt(2);
	U[1][1] = -1 / sqrt(2);
	omp_set_num_threads(n_threads);
	double timer1 = omp_get_wtime();
	complexd * a = generate(n);
	timer1 = omp_get_wtime() - timer1;
	double timer2 = omp_get_wtime();
	complexd * b = quantum(a, n, U, k);
	timer2 = omp_get_wtime() - timer2;
	long long size = power(2, n);
	ofstream f("vectors.txt");
	for (long long i = 0; i < size; i++) {
		f << a[i] << "\n";
	}
	f << endl;
	for (long long i = 0; i < size; i++) {
		f << b[i] << "\n";
	}
	//cout << timer1 << endl;
	cout << timer2 << endl;
	//cin >> n;
	return 0;
}