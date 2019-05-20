TARGET=task5 print
all: $(TARGET)

#TODO add entries for different tests

task5: task5.cpp quantum.h
	mpicxx -O3 -Wall -Werror -std=c++11 task5.cpp -o task5 -fopenmp
print: print.cpp
	g++ -O3 print.cpp -o print
clean:
	rm -rf $(TARGET) *.o
test: task5 print
	mpirun -n 2 ./task5 2 b.dat a_test_1.dat test
	diff -s b.txt b_test_1.txt

	mpirun -n 2 ./task5 2 b.dat a_test_2.dat test
	diff -s b.txt b_test_2.txt

	mpirun -n 2 ./task5 3 b.dat a_test_3.dat test
	diff -s b.txt b_test_3.txt

lint: task5
	cpplint --filter=-legal/copyright,-build/namespaces task5.cpp
	cpplint --filter=-whitespace,-legal/copyright,-build/namespaces --linelength=120 quantum.h
