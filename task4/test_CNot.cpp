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
    complexd *source = new complexd[size];
    complexd *transformed = new complexd[size];
    if (fread(source, sizeof(source[0]), size, f) != size) {
        //exit(1);
    }
    fclose(f);

    int n = atoi(argv[2]);
    k = n - k;
    l = n - l;
    uint64_t bitk = 1ull << k;
    uint64_t bitl = 1ull << l;
    for (uint64_t i = 0; i < size; ++i) {
        if (i & bitk) {
            transformed[i] = source[i ^ bitl];
        } else {
            transformed[i] = source[i];
        }
    }

    f = fopen(argv[4], "wb");
    fwrite(&size, sizeof(size), 1, f);
    fwrite(transformed, sizeof(transformed[0]), size, f);
    fclose(f);

    delete [] source;
    delete [] transformed;
}
