.PHONY: all clean

CXX=g++
CXXFLAGS=-O3 -std=c++11 -fopenmp -g

xg=../bin/xg

all: test clean

test: $(xg)
	prove -v t

$(xg):
	cd .. && $(MAKE) bin/xg

#clean:
#	rm -f build_graph
