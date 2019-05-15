TARGET=task5 print
all: $(TARGET)

#TODO add entries for different tests

task5: task5.cpp quantum.h
	mpicxx -O3 -Wall -Werror -std=c++11 task5.cpp -o task5 -fopenmp
print: print.cpp
	g++ -O3 print.cpp -o print
clean:
	rm -rf $(TARGET) *.o
#
#test: task4 print Not Adamar CNot
#	mpirun -n 2 ./task4 Not 1 1 out.txt a_out.txt
#	./test_not a_out.txt 1 1 test_out.txt
#	./print test_out.txt 1 > test_out
#	./print out.txt 1 > out
#	diff -s out test_out
#
#	mpirun -n 2 ./task4 Adamar 3 1 out.txt a_out.txt
#	./test_Adamar out.txt 3 1 test_out.txt
#	./print test_out.txt 1 > test_out
#	./print a_out.txt 1 > out
#	diff -s out test_out
#
#	mpirun -n 2 ./task4 CNot 2 1 2 out.txt a_out.txt
#	./test_CNot a_out.txt 2 1 2 test_out.txt
#	./print test_out.txt 2 > test_out
#	./print out.txt 2 > out
#	diff -s out test_out
#
#	# w=0
#	mpirun -n 2 ./task4 Rotate 1 2 0 out.txt a_out.txt
#	./test_not a_out.txt 1 2 test_out.txt
#	./print test_out.txt 1 > test_out
#	./print out.txt 1 > out
#	diff -s out test_out
#
#	mpirun -n 2 ./task4 CRotate 2 1 2 "pi/2" out.txt a_out.txt
#	./print a_out.txt 2 > test_out
#	./print out.txt 2 > out
#	diff -s out test_out
#

lint: task5
	cpplint --filter=-legal/copyright,-build/namespaces task5.cpp
	cpplint --filter=-whitespace,-legal/copyright,-build/namespaces --linelength=120 quantum.h
