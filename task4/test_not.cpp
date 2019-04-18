#include <iostream>
#include <complex>

void usage(int argc, char **argv) {
    printf("Usage: %s <input_file> <k> <output_file>\n"
            "   k - qubit index\n"
            "File format:\n"
            "   uint64_t n, //vector size, not number of qubits\n"
            "   complex<double> data[n]\n", argv[0]);
    exit(1);
}

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

int main(int argc, char **argv) {
    // if (argc != 4) {
    //     usage(argc, argv);
    // }
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
    std::complex<double> *a = new std::complex<double>[size];
    std::complex<double> *b = new std::complex<double>[size];
    for (long long i = 0; i < size; i++) {
        if (fread(&a[i], sizeof(a[0]), 1, f) != 1) {
            //exit(1);
        }
    }
    fclose(f);

    int q = -1;
    uint64_t tmpn = size;
    while (tmpn) {
        q++;
        tmpn >>= 1;
    }
    k = q - k;
    uint64_t bit = 1ull << k;
    for (uint64_t i = 0; i < size; ++i) {
        b[i] = a[i ^ bit];
    }

    f = fopen(argv[4], "wb");
    fwrite(&size, sizeof(size), 1, f);
    fwrite(&k, sizeof(k), 1, f);
    fwrite(b, sizeof(b[0]), size, f);
    fclose(f);

    delete [] a;
    delete [] b;
}
