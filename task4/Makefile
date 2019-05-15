TARGET=task4 print
all: $(TARGET)

#TODO add entries for different tests

task4: task4.cpp quantum.h
	mpicxx -O3 -Wall -Werror -std=c++11 task4.cpp -o task4 -fopenmp
print: print.cpp
	g++ -O3 print.cpp -o print
clean:
	rm -rf $(TARGET) *.o
Not: test_not.cpp
	g++ test_not.cpp -o test_not
Adamar: test_Adamar.cpp
	g++ test_Adamar.cpp -o test_Adamar
CNot: test_CNot.cpp
	g++ test_CNot.cpp -o test_CNot

test: task4 print Not Adamar CNot
	mpirun -n 2 ./task4 Not 1 1 out.txt a_out.txt
	./test_not a_out.txt 1 1 test_out.txt
	./print test_out.txt 1 > test_out
	./print out.txt 1 > out
	diff -s out test_out

	mpirun -n 2 ./task4 Adamar 3 1 out.txt a_out.txt
	./test_Adamar out.txt 3 1 test_out.txt
	./print test_out.txt 1 > test_out
	./print a_out.txt 1 > out
	diff -s out test_out

	mpirun -n 2 ./task4 CNot 2 1 2 out.txt a_out.txt
	./test_CNot a_out.txt 2 1 2 test_out.txt
	./print test_out.txt 2 > test_out
	./print out.txt 2 > out
	diff -s out test_out

	# w=0
	mpirun -n 2 ./task4 Rotate 1 2 0 out.txt a_out.txt
	./test_not a_out.txt 1 2 test_out.txt
	./print test_out.txt 1 > test_out
	./print out.txt 1 > out
	diff -s out test_out

	mpirun -n 2 ./task4 CRotate 2 1 2 "pi/2" out.txt a_out.txt
	./print a_out.txt 2 > test_out
	./print out.txt 2 > out
	diff -s out test_out


lint: task4
	cpplint --filter=-legal/copyright,-build/namespaces task4.cpp 
	cpplint --filter=-legal/copyright,-build/namespaces test_not.cpp 
	cpplint --filter=-legal/copyright,-build/namespaces test_Adamar.cpp
	cpplint --filter=-whitespace,-legal/copyright,-build/namespaces --linelength=120 quantum.h
