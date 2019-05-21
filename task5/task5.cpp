#include <string>
#include <fstream>
#include "./quantum.h"

uint64_t invert_index(uint64_t idx, int n) {
    uint64_t bit, inverted = 0;
    for (int i = 0; i < n; i++)
    {
        bit = (idx >> i) & 1;
        inverted = (inverted << 1) | bit;
    }
    //cout << inverted << endl;
    return inverted;
}

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);
    //omp_set_num_threads(1);
    int n = atoi(argv[1]);
    bool test = false;
    if (argc == 5) test = true;

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    uint64_t size_array = power(2, n);
    complexd * a, * a_file;
    if (!test) a = generate(n);
    if (test) {
        if (rank == 0) cout << "test" << endl;
        MPI_File fh;
        int complex_size, int_size;
        MPI_Type_size(MPI_CXX_DOUBLE_COMPLEX, &complex_size);
        MPI_Type_size(MPI_LONG_LONG_INT, &int_size);
        MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        MPI_File_read_all(fh, &size_array, 1, MPI_LONG_LONG_INT, MPI_STATUS_IGNORE);
        MPI_File_set_view(fh, rank * size_array/size * complex_size + int_size, MPI_CXX_DOUBLE_COMPLEX,
                          MPI_CXX_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
        a_file = new complexd[size_array/size];
        MPI_File_read(fh, a_file, size_array/size, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        normalize(a_file, size_array/size);
        MPI_File_close(&fh);
        a = new complexd[size_array/size];
        memmove(a, a_file, sizeof(complexd) * size_array/size);
    }
    File_MPI_Write(a, size_array, 0, 0, argv[3], rank);
    complexd * b = new complexd[size_array/size];
    double timer = MPI_Wtime();
    for (int i = 1; i <= n; i++) {
        //cout << i << endl;
        //double timer2 = MPI_Wtime();
        b = quantum_Adamar(a, n, i);
        //timer2 = MPI_Wtime() - timer2;
        //timer += timer2;
        memmove(a, b, sizeof(complexd) * size_array / size);
        int m = 2;
        for (int j = i + 1; j <= n; j++) {
            //if (rank == 0) cout << i << " " << j << endl;
            double phi = 2 * M_PI / power(2, m);
            //timer2 = MPI_Wtime();
            b = quantum_CRotate(a, n, i, j, phi);
            //timer2 = MPI_Wtime() - timer2;
            //timer += timer2;
            memmove(a, b, sizeof(complexd) * size_array / size);
            m++;
            if (j != n) delete[] b;
        }
        if (i != n) delete[] b;
        //cout << endl;
    }
    timer = MPI_Wtime() - timer;
    if (rank == 0) cout << timer << endl;
    File_MPI_Write(b, size_array, 0, 0, argv[2], rank);
    if (test) {
        complexd * b_out = new complexd[size_array];
        MPI_Gather(b, size_array/size, MPI_CXX_DOUBLE_COMPLEX, b_out, size_array/size, MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            fstream f("b.txt", ios::out);
            f << size_array << endl;
            for (uint64_t i = 0; i < size_array; i++) {
                uint64_t idx = invert_index(i, n);
                f << b_out[idx] << endl;
            }
            f.close();
        }
    }
    delete[] a;
    delete[] b;
    MPI_Finalize();
    return 0;
}
