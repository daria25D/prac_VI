#include <string>
#include <fstream>
#include "./quantum.h"

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);
    omp_set_num_threads(1);
    int n = atoi(argv[1]);
    complexd * a = generate(n);
    uint64_t size_array = power(2, n);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //File_MPI_Write(a, size_array, 0, 0, argv[3], rank);
    complexd * b = new complexd[size_array/size];
    double timer = MPI_Wtime();
    for (int i = 1; i <= n; i++) {
        //cout << i << endl;
        double timer2 = MPI_Wtime();
        b = quantum_Adamar(a, n, i);
        timer2 = MPI_Wtime() - timer2;
        timer += timer2;
        memmove(a, b, sizeof(complexd) * size_array / size);
        int m = 2;
        for (int j = i + 1; j <= n; j++) {
            //if (rank == 0) cout << i << " " << j << endl;
            double phi = 2 * M_PI / power(2, m);
            timer2 = MPI_Wtime();
            b = quantum_CRotate(a, n, i, j, phi);
            timer2 = MPI_Wtime() - timer2;
            timer += timer2;
            memmove(a, b, sizeof(complexd) * size_array / size);
            m++;
            if (j != n) delete[] b;
        }
        if (i != n) delete[] b;
        //cout << endl;
    }
    timer = MPI_Wtime() - timer;
    if (rank == 0) cout << timer << endl;
    //File_MPI_Write(b, size_array, 0, 0, argv[2], rank);
    //File_MPI_Write(a, size_array, 0, 0, argv[3], rank);
    delete[] a;
    delete[] b;
    MPI_Finalize();
    return 0;
}
