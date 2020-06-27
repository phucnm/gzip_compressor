EXTRA_CXXFLAGS= -g
EXTRA_CFLAGS=
CXXFLAGS=-O3 -Wall -std=c++17 $(EXTRA_CXXFLAGS)
CFLAGS=-O3 -Wall -std=c11 $(EXTRA_CFLAGS)

all: uvgz

clean:
	rm -f uvgz *.o
