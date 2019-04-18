#include <iostream>
#include <complex>

typedef std::complex<double> complexd;
using namespace std;

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

int main(int argc, char **argv) {
    int n, k;
    n = atoi(argv[2]);
    k = atoi(argv[3]);
    uint64_t size;
    FILE *f = fopen(argv[1], "rb");
    if (fread(&size, sizeof(size), 1, f) != 1) {
        cout << 1 << endl;
        exit(1);
    }
    if (fread(&k, sizeof(k), 1, f) != 1) {
        cout << 2 << endl;
        exit(1);
    }
    complexd *a = new complexd[size];
    complexd *b = new complexd[size];
    for (uint64_t i = 0; i < size; i++) {
        if (fread(&a[i], sizeof(a[0]), 1, f) != 1) {
            // exit
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
