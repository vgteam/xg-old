.PHONY: all clean

CXX=g++
CXXFLAGS=-O3 -std=c++11 -fopenmp -g

xg=../xg

all: test clean

test: $(xg)
	prove -v t

$(xg):
	cd .. && $(MAKE) xg

#clean:
#	rm -f build_graph
