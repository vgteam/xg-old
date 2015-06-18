.PHONY: all clean

CXX=g++
CXXFLAGS=-O3 -std=c++11 -fopenmp -g

succinctg=../succinctg

all: test clean

test: $(succinctg)
	prove -v t

$(succinctg):
	cd .. && $(MAKE) succinctg

#clean:
#	rm -f build_graph
