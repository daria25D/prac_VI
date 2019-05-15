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
    // TODO redo for different cases (like CNot, CRotate)
    //./print what > where
    if (argc != 3) {
        cerr << "Wrong number of arguments";
        return -1;
    }
    fstream file(argv[1], ios::binary | ios::in);
    int type = atoi(argv[2]);
    long long size;
    int k, l;
    file.read((char *)&size, sizeof(size));
    file.read((char *)&k, sizeof(k));
    if (type == 2) {
        file.read((char *)&l, sizeof(l));
    }
    cout << size << setw(5) << k;
    if (type == 2) {
        cout << setw(5) << l;
    }
    cout << endl;
    complex<double> f;
    for (long long i = 0; i < size; i++) {
        file.read((char *)&f, sizeof(f));
        cout << f << endl;
        
    }
    return 0;
}