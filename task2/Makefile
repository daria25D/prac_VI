TARGET=task2 gen print
all: $(TARGET)

task2: task2.cpp
	mpicxx -O3 -std=c++11 task2.cpp -o task2
gen: gen.cpp
	g++ -O3 gen.cpp -o gen
print: print.cpp
	g++ -O3 print.cpp -o print
clean:
	rm -rf $(TARGET) *.o