#include "quantum.h"
#include <cstring>
#include <fstream>

int main(int argc, char ** argv) {
	string type = string(argv[1]);
		MPI_Init(&argc, &argv);
	//fstream f("out.txt", ios::out);
	if (type == "Adamar") {
		//Adamar n k
		int n = atoi(argv[2]);
		int k = atoi(argv[3]);
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		complexd * a = generate(n);
		complexd * b = quantum_Adamar(a, n, k);
		File_MPI_Write(b, power(2, n), k, argv[4], rank);
		delete[] a;
		delete[] b;
	} else if (type == "nAdamar") {
		//nAdamar n
		int n = atoi(argv[2]);
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		complexd * a = generate(n);
		complexd * b = quantum_nAdamar(a, n);
		File_MPI_Write(b, power(2, n), 0, argv[3], rank);
		delete[] a;
		delete[] b;

	} else if (type == "Rotate") {
		//Rotate n k phi
		int n = atoi(argv[2]);
		int k = atoi(argv[3]);
		double phi = atof(argv[4]);
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		complexd * a = generate(n);
		complexd * b = quantum_Rotate(a, n, k, phi);
		File_MPI_Write(b, power(2, n), k, argv[5], rank);
		delete[] a;
		delete[] b;

	} else if (type == "Not") {
		//Not n k
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

	} else {

	}
	MPI_Finalize();
	//f.close();
	return 0;
}