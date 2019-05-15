#include <iostream>
#include <complex>


typedef std::complex<double> complexd;
using namespace std;


int main(int argc, char **argv) {

    int k, l;
    sscanf(argv[3], "%d", &k);
    sscanf(argv[4], "%d", &l);
    uint64_t size;
    FILE *f = fopen(argv[1], "rb");
    if (fread(&size, sizeof(size), 1, f) != 1) {
        exit(1);
    }
    if (fread(&k, sizeof(k), 1, f) != 1) {
        cout << 2 << endl;
        exit(1);
    }
    if (fread(&l, sizeof(l), 1, f) != 1) {
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
    int q = k, j = l;
    int n = atoi(argv[2]);
    k = n - k;
    l = n - l;
    uint64_t bitk = 1ull << k;
    uint64_t bitl = 1ull << l;
    for (uint64_t i = 0; i < size; i++) {
        if (i & bitk) {
            b[i] = a[i ^ bitl];
        } else {
            b[i] = a[i];
        }
    }
//    for (uint64_t i = 0; i < size; i++) {
//        cout << b[i] << endl;
//    }
    f = fopen(argv[5], "wb");
    fwrite(&size, sizeof(size), 1, f);
    fwrite(&q, sizeof(q), 1, f);
    fwrite(&j, sizeof(j), 1, f);
    fwrite(b, sizeof(b[0]), size, f);
    fclose(f);

    delete [] a;
    delete [] b;
}
