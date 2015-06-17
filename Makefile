.PHONY: all clean test get-deps

CXX=g++
CXXFLAGS=-O3 -std=c++11 -fopenmp -g
LIBS=main.o cpp/vg.pb.o
INCLUDES=-I./ -Icpp -Istream/protobuf/build/include -Isdsl-lite/build/include -Istream
LDFLAGS=-L./ -Lstream/protobuf -Lsdsl-lite/build/lib -lprotobuf -lsdsl -lz
STREAM=stream
PROTOBUF=$(STREAM)/protobuf
LIBPROTOBUF=stream/protobuf/libprotobuf.a
LIBSDSL=sdsl-lite/build/lib/libsdsl.a
EXECUTABLE=succinctg

all: $(EXECUTABLE) doc

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

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp $(INCLUDES)

$(EXECUTABLE): $(LIBS)
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $(LIBS) $(INCLUDES) $(LDFLAGS)

clean:
	rm -rf cpp
	rm -f $(EXECUTABLE)
	rm -f *.o
	rm $(LIBPROTOBUF)
	cd $(PROTOBUF) && $(MAKE) clean && rm -rf build
