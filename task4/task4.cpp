/* "Copyright 2019 <Copyright Owner>" */
#include <string>
#include <fstream>
#include "./quantum.h"

int main(int argc, char ** argv) {
    string type = string(argv[1]);
    MPI_Init(&argc, &argv);
    if (type == "Adamar") {
        // Adamar n k file_b file_a
        int n = atoi(argv[2]);
        int k = atoi(argv[3]);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        complexd * a = generate(n);
        complexd * b = quantum_Adamar(a, n, k);
        File_MPI_Write(b, power(2, n), k, argv[4], rank);
        File_MPI_Write(a, power(2, n), k, argv[5], rank);
        delete[] a;
        delete[] b;
    } else if (type == "nAdamar") {
        // nAdamar n file_b file_a
        int n = atoi(argv[2]);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        complexd * a = generate(n);
        complexd * b = quantum_nAdamar(a, n);
        File_MPI_Write(b, power(2, n), 0, argv[3], rank);
        File_MPI_Write(a, power(2, n), 0, argv[4], rank);
        delete[] a;
        delete[] b;
    } else if (type == "Rotate") {
        // Rotate n k phi file_b file_a
        int n = atoi(argv[2]);
        int k = atoi(argv[3]);
        double phi = atof(argv[4]);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        complexd * a = generate(n);
        complexd * b = quantum_Rotate(a, n, k, phi);
        File_MPI_Write(b, power(2, n), k, argv[5], rank);
        File_MPI_Write(a, power(2, n), k, argv[6], rank);
        delete[] a;
        delete[] b;
    } else if (type == "Not") {
        // Not n k file_b file_a
        int n = atoi(argv[2]);
        int k = atoi(argv[3]);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        complexd * a = generate(n);
        complexd * b = quantum_Not(a, n, k);
        File_MPI_Write(b, power(2, n), k, argv[4], rank);
        File_MPI_Write(a, power(2, n), k, argv[5], rank);
        delete[] a;
        delete[] b;
    } else if (type == "CNot") {
        // CNot n k l file_b file_a
        int n = atoi(argv[2]);
        int k = atoi(argv[3]);
        int l = atoi(argv[4]);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        complexd * a = generate(n);
        complexd * b = quantum_CNot(a, n, k, l);
        File_MPI_Write(b, power(2, n), k, argv[5], rank);
        File_MPI_Write(a, power(2, n), k, argv[6], rank);
        delete[] a;
        delete[] b;
    } else if (type == "CRotate") {
        // CRotate n k l phi file_b file_a
        int n = atoi(argv[2]);
        int k = atoi(argv[3]);
        int l = atoi(argv[4]);
        double phi = atof(argv[5]);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        complexd * a = generate(n);
        complexd * b = quantum_CRotate(a, n, k, l, phi);
        File_MPI_Write(b, power(2, n), k, argv[6], rank);
        File_MPI_Write(a, power(2, n), k, argv[7], rank);
        delete[] a;
        delete[] b;
    }
    MPI_Finalize();
    return 0;
}
