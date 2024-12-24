CC=gcc
CXX=g++
CFLAGS=-Wall -Wextra -g -fopenmp -Os
CXXFLAGS=-Wall -Wextra -g -fopenmp -Os
LDFLAGS=-lm

all : dos dos_cpp test

dos: dos.c
	$(CC) $(CFLAGS) dos.c -o dos $(LDFLAGS)

dos_cpp: dos.cpp
	$(CXX) $(CXXFLAGS) dos.cpp -o dos_cpp $(LDFLAGS)

test: test.c
	$(CC) $(CFLAGS) test.c -o test $(LDFLAGS)

clean:
	$(RM) dos dos_cpp test
