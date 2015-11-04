.PHONY: all clean test

LIB_DIR:=lib
OBJ_DIR:=obj
BIN_DIR:=bin
SRC_DIR:=src
CPP_DIR:=cpp

CXX=g++
CXXFLAGS=-O3 -std=c++11 -fopenmp -g
OBJ=cpp/vg.pb.o xg.o # main.o not included for easier libxg.a creation
LD_INCLUDES=-I./ -Icpp -Istream -I$(SRC_DIR)
LD_LIBS=-lprotobuf -lsdsl -lz -ldivsufsort -ldivsufsort64 -lgomp -lm -lpthread
STREAM=stream
EXE:=xg


#Some little adjustments to build on OSX
#(tested with gcc4.9 installed from MacPorts)
SYS=$(shell uname -s)
all: $(EXE)

doc: README.md
README.md: README.base.md
	pandoc -o README.html -s README.base.md
	pandoc -o DESIGN.html -s DESIGN.md
	cat README.base.md >README.md
	cat DESIGN.html | tail -7| perl -p -e 's/<p>/\n/g' | sed 's%</p>%%g' | head -10 >>README.md



$(CPP_DIR)/vg.pb.cc: $(CPP_DIR)/vg.pb.h
$(CPP_DIR)/vg.pb.h: $(SRC_DIR)/vg.proto
	mkdir -p cpp
	protoc $(SRC_DIR)/vg.proto --cpp_out=cpp

$(OBJ_DIR)/vg.pb.o: $(CPP_DIR)/vg.pb.h $(CPP_DIR)/vg.pb.cc
	$(CXX) $(CXXFLAGS) -c -o cpp/vg.pb.o cpp/vg.pb.cc $(INCLUDES)

$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp $(CPP_DIR)/vg.pb.h $(SRC_DIR)/xg.hpp 
	$(CXX) $(CXXFLAGS) $(LD_LIBS) -c -o $(OBJ_DIR)/main.o $(SRC_DIR)/main.cpp $(LD_INCLUDES)

$(OBJ_DIR)/xg.o: $(SRC_DIR)/xg.cpp $(SRC_DIR)/xg.hpp $(CPP_DIR)/vg.pb.h
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDES)

$(EXE): $(OBJ_DIR)/main.o
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $< $(LD_INCLUDES) $(LD_LIBS) $(STATICFLAGS)

libxg.a: $(LIBS)
	ar rs libxg.a $(LIBS)

test:
	cd test && make

clean-xg:
	rm -f *.o xg cpp/* libxg.a

clean:
	rm -rf cpp
	rm -f $(EXECUTABLE)
	rm -f *.o
	rm -f libxg.a
