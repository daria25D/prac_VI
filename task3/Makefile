TARGET=task3 gen print
all: $(TARGET)

task3: task3.cpp
	mpicxx -O3 -std=c++11 task3.cpp -o task3 -fopenmp
gen: gen.cpp
	g++ -O3 gen.cpp -o gen
print: print.cpp
	g++ -O3 print.cpp -o print
clean:
	rm -rf $(TARGET) *.o