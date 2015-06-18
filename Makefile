.PHONY: all clean test get-deps

CXX=g++
CXXFLAGS=-O3 -std=c++11 -fopenmp -g
LIBS=main.o cpp/vg.pb.o succinct.o
INCLUDES=-I./ -Icpp -Istream/protobuf/build/include -Isdsl-lite/build/include -Isdsl-lite/build/external/libdivsufsort-2.0.1/include -Istream
LDFLAGS=-L./ -Lstream/protobuf -Lsdsl-lite/build/lib -lprotobuf -lsdsl -lz
STREAM=stream
PROTOBUF=$(STREAM)/protobuf
LIBPROTOBUF=stream/protobuf/libprotobuf.a
LIBSDSL=sdsl-lite/build/lib/libsdsl.a
EXECUTABLE=succinctg

all: $(EXECUTABLE)

doc: README.md
	pandoc -o README.html -s README.md
	pandoc -o README.pdf README.md

$(LIBPROTOBUF):
	cd $(STREAM) && $(MAKE)

$(LIBSDSL):
	cd sdsl-lite/build && cmake .. -Wno-dev && $(MAKE)

cpp/vg.pb.cc: cpp/vg.pb.h
cpp/vg.pb.h: vg.proto $(LIBPROTOBUF)
	mkdir -p cpp
	$(PROTOBUF)/build/bin/protoc vg.proto --cpp_out=cpp
cpp/vg.pb.o: cpp/vg.pb.h cpp/vg.pb.cc
	$(CXX) $(CXXFLAGS) -c -o cpp/vg.pb.o cpp/vg.pb.cc $(INCLUDES)

main.o: main.cpp $(LIBSDSL) cpp/vg.pb.h
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp $(INCLUDES)

succinct.o: succinct.cpp succinct.hpp $(LIBSDSL) cpp/vg.pb.h
	$(CXX) $(CXXFLAGS) -c -o succinct.o succinct.cpp $(INCLUDES)

$(EXECUTABLE): $(LIBS)
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $(LIBS) $(INCLUDES) $(LDFLAGS)

get-deps:
	echo todo

test:
	cd test && make

clean:
	rm -rf cpp
	rm -f $(EXECUTABLE)
	rm -f *.o
	rm $(LIBPROTOBUF)
	cd $(PROTOBUF) && $(MAKE) clean && rm -rf build
