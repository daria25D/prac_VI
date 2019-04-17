#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>

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

int main(int argc, char ** argv) {
    //./print what > where
    if (argc != 2) {
        cerr << "Wrong number of arguments";
        return -1;
    }
    fstream file(argv[1], ios::binary | ios::in);
    int n, k;
    file.read((char *)&n, sizeof(n));
    file.read((char *)&k, sizeof(k));
    cout << n << setw(5) << k << endl;
    complex<double> f;
    long long size = power(2, n);
    for (long long i = 0; i < size; i++) {
        file.read((char *)&f, sizeof(f));
        cout << f << endl;
        
    }
    return 0;
}