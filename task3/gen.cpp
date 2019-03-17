#include <iostream>
#include <fstream>
#include <complex>
#include <random>

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

int main(int argc, char **argv) {
    //./generate n file
    if (argc != 3) {
        cerr << "Wrong number of arguments";
        return -1;
    }
    srand(time(NULL));
    int n;
    try {
        n = strtoull(argv[1], NULL, 0);
        //k = strtoull(argv[2], NULL, 0);
    } catch (exception & e) {
        cerr << e.what();
        return -1;
    }
    fstream file(argv[2], ios::binary | ios::out);
    if (!file.is_open()) {
        cerr << "Cannot open file";
        return -1;
    }
    file.write((char *)&n, sizeof(n));
    //file.write((char *)&k, sizeof(k));
	unsigned int seed = time(0);
    long long size = pow(2, n);
    try {
        for (int i = 0; i < size; i++) {
            complex<double> f;
            f.real(rand_r(&seed));
            f.imag(rand_r(&seed));
            file.write((char *)&f, sizeof(complex<double>));
        }            
    } catch (exception & e) {
        cerr << e.what();
        return -1;
    } 
    return 0;
} 