.PHONY: all clean test get-deps

CXX=g++
CXXFLAGS=-O3 -std=c++11 -fopenmp -g
LIBS=cpp/vg.pb.o xg.o # main.o not included for easier libxg.a creation
INCLUDES=-I./ -Icpp -Istream/protobuf/build/include -Isdsl-lite/build/include -Isdsl-lite/build/external/libdivsufsort-2.0.1/include -Istream
LDSEARCH=-L./ -Lstream/protobuf -Lsdsl-lite/build/lib -Lsdsl-lite/build/external/libdivsufsort-2.0.1/lib
LDFLAGS=-lprotobuf -lsdsl -lz -ldivsufsort -ldivsufsort64 -lgomp -lm -lpthread
STREAM=stream
PROTOBUF=$(STREAM)/protobuf
LIBPROTOBUF=stream/protobuf/libprotobuf.a
LIBSDSL=sdsl-lite/build/lib/libsdsl.a
EXECUTABLE=xg

#Some little adjustments to build on OSX
#(tested with gcc4.9 installed from MacPorts)
SYS=$(shell uname -s)
ifeq (${SYS},Darwin)
	CMAKE_BIN= # linux binary won't work on os x
	CMAKE_SETPATH= # mac ports cmake seems to work just fine so leave blank
	STATICFLAGS= # -static doesn't work on OSX unless libgcc compiled as static.
else
	CMAKE_BIN=cmake-3.3.0-rc2-Linux-x86_64/bin/cmake
	CMAKE_SETPATH=PATH=../../cmake-3.3.0-rc2-Linux-x86_64/bin/:${PATH}
	STATICFLAGS=-static -static-libstdc++ -static-libgcc -Wl,-Bstatic
endif

all: $(EXECUTABLE)

doc: README.md
README.md: README.base.md
	pandoc -o README.html -s README.base.md
	pandoc -o DESIGN.html -s DESIGN.md
	cat README.base.md >README.md
	cat DESIGN.html | tail -7| perl -p -e 's/<p>/\n/g' | sed 's%</p>%%g' | head -10 >>README.md

$(LIBPROTOBUF):
	cd $(STREAM) && $(MAKE)

$(LIBSDSL): $(CMAKE_BIN)
	$(CMAKE_SETPATH) cmake --version
	cd sdsl-lite/build && $(CMAKE_SETPATH) cmake .. -Wno-dev && $(MAKE)

cpp/vg.pb.cc: cpp/vg.pb.h
cpp/vg.pb.h: vg.proto $(LIBPROTOBUF)
	mkdir -p cpp
	$(PROTOBUF)/build/bin/protoc vg.proto --cpp_out=cpp
cpp/vg.pb.o: cpp/vg.pb.h cpp/vg.pb.cc
	$(CXX) $(CXXFLAGS) -c -o cpp/vg.pb.o cpp/vg.pb.cc $(INCLUDES)

main.o: main.cpp $(LIBSDSL) cpp/vg.pb.h xg.hpp 
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp $(INCLUDES)

xg.o: xg.cpp xg.hpp $(LIBSDSL) cpp/vg.pb.h
	$(CXX) $(CXXFLAGS) -c -o xg.o xg.cpp $(INCLUDES)

$(EXECUTABLE): $(LIBS) main.o
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $(LIBS) main.o $(INCLUDES) $(LDSEARCH) $(STATICFLAGS) $(LDFLAGS)

libxg.a: $(LIBS)
	ar rs libxg.a $(LIBS)

$(CMAKE_BIN):
	wget http://www.cmake.org/files/v3.3/cmake-3.3.0-rc2-Linux-x86_64.tar.gz
	tar xzvf cmake-3.3.0-rc2-Linux-x86_64.tar.gz

get-deps: $(CMAKE_BIN)


test:
	cd test && make

clean-xg:
	rm -f *.o xg cpp/* libxg.a

clean:
	rm -rf cpp
	rm -f $(EXECUTABLE)
	rm -f *.o
	rm -f libxg.a
	rm $(LIBPROTOBUF)
	cd $(PROTOBUF) && $(MAKE) clean && rm -rf build
	rm -rf cmake-3.3.0-rc2-Linux-x86_64.tar.gz cmake-3.3.0-rc2-Linux-x86_64
