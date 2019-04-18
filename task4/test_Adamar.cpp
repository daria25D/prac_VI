#include <iostream>
#include <complex>

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

std::complex<double> * quantum(std::complex<double> * a, long long n, std::complex<double> U[][2], long long k) {
	long long size_of_a = power(2, n);
	std::complex<double> * b = new std::complex<double>[size_of_a];
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

int main(int argc, char **argv) {
    int n, k;
    n = atoi(argv[2]);
    k = atoi(argv[3]);
    long long size;
    FILE *f = fopen(argv[1], "rb");
    if (fread(&size, sizeof(size), 1, f) != 1) {
        std::cout << 1 << std::endl;
        exit(1);
    }
    if (fread(&k, sizeof(k), 1, f) != 1) {
        std::cout << 2 << std::endl;
        exit(1);
    }
    std::complex<double> *b = new std::complex<double>[size];
    std::complex<double> *a = new std::complex<double>[size];
    for (long long i = 0; i < size; i++) {
        if (fread(&b[i], sizeof(b[0]), 1, f) != 1) {
            //exit(1);
        }
    }
    fclose(f);
    std::complex<double> U[2][2] = {
                            {1 / sqrt(2), 1 / sqrt(2)}, 
                            {1 / sqrt(2), -1 / sqrt(2)}};

    a = quantum(b, n, U, k);
    f = fopen(argv[4], "wb");
    fwrite(&size, sizeof(size), 1, f);
    fwrite(&k, sizeof(k), 1, f);
    fwrite(a, sizeof(a[0]), size, f);
    fclose(f);

    delete [] b;
    delete [] a;
}
