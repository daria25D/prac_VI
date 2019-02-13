TARGET=task1
all: $(TARGET)

task1: task1.cpp
	g++ -O0 -std=c++11 task1.cpp -o task1 -fopenmp

clean:
	rm -rf $(TARGET) *.o