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

complexd * quantum(complexd * a, uint64_t n, complexd U[][2], uint64_t k) {
    uint64_t size_of_a = power(2, n);
    complexd * b = new complexd[size_of_a];
    for (uint64_t i = 0; i < size_of_a; i++) {
        uint64_t idx1, idx2;
        uint64_t idx_ik = (i >> (n - k)) % 2;
        uint64_t n1 = 1;
        for (uint64_t j = 0; j < n - k; j++)
            n1 <<= 1;
        uint64_t n0 = 1;
        for (uint64_t j = 0; j < n; j++) {
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
    complexd *b = new complexd[size];
    complexd *a = new complexd[size];
    for (uint64_t i = 0; i < size; i++) {
        if (fread(&b[i], sizeof(b[0]), 1, f) != 1) {
            // exit(1);
        }
    }
    fclose(f);
    complexd U[2][2] = {
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
