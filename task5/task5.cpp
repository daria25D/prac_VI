#include <string>
#include <fstream>
#include "./quantum.h"

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);
    int n = atoi(argv[1]);
    complexd * a = generate(n);
    uint64_t size_array = power(2, n);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    complexd * b = new complexd[size_array/size];
    for (uint64_t i = 0; i < size_array; i++) {
        b = quantum_Adamar(a, n, i + 1);
        int m = 2;
        for (uint64_t j = i + 1; j < size_array; j++) {
            double phi = 2 * M_PI / power(2, m);
            b = quantum_CRotate(b, n, i + 1, j, phi);
            m++;
        }
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    File_MPI_Write(b, size_array, 0, 0, argv[2], rank);
    File_MPI_Write(a, size_array, 0, 0, argv[3], rank);
    delete[] a;
    delete[] b;
    MPI_Finalize();
    return 0;
}
