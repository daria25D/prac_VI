TARGET=task4 print
all: $(TARGET)

task4: task4.cpp quantum.h
	mpicxx -O3 -Wall -Werror -std=c++11 task4.cpp -o task4 -fopenmp
print: print.cpp
	g++ -O3 print.cpp -o print
clean:
	rm -rf $(TARGET) *.o

test: task4
	mpirun -n 2 ./task4 Not 2 1 out.txt a_out.txt
	g++ test_not.cpp -o test_not
	./test_not a_out.txt 2 1 test_out.txt
	diff -s out.txt test_out.txt

lint: task4
	.cpplint task4.cpp 
