CXX = g++-11 # Compiler to use
CFLAGS = -fopenmp -g -fdiagnostics-color=always # Flags to pass


output: main.o
	$(CXX) main.o -o main $(CFLAGS)

main.o: main.cpp
	$(CXX) -c main.cpp

clean:
	rm *.o main