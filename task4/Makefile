TARGET=task4 print
all: $(TARGET)

task4: task4.cpp quantum.h
	mpicxx -O3 -Wall -Werror -std=c++11 task4.cpp -o task4 -fopenmp
print: print.cpp
	g++ -O3 print.cpp -o print
clean:
	rm -rf $(TARGET) *.o

test: task4 print
	mpirun -n 2 ./task4 Not 2 1 out.txt a_out.txt
	g++ test_not.cpp -o test_not
	./test_not a_out.txt 2 1 test_out.txt
	./print test_out.txt > test_out
	./print out.txt > out
	diff -s out test_out
	mpirun -n 2 ./task4 Adamar 3 1 out.txt a_out.txt
	g++ test_Adamar.cpp -o test_Adamar
	./test_Adamar out.txt 3 1 test_out.txt
	./print test_out.txt > test_out
	./print a_out.txt > out
	diff -s out test_out

lint: task4
	cpplint --filter=-legal/copyright,-build/namespaces task4.cpp 
	cpplint --filter=-legal/copyright,-build/namespaces test_not.cpp 
	cpplint --filter=-legal/copyright,-build/namespaces test_Adamar.cpp 
