MouseDNA: MouseDNA.o
	g++ -Wall -O -g -std=c++11 -o MouseDNA MouseDNA.o -fopenmp
MouseDNA.o: MouseDNA.h MouseDNA.cpp
	g++ -Wall -O -g -std=c++11 -c -o MouseDNA.o MouseDNA.cpp -fopenmp
clean:
	rm *.o
