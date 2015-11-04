.PHONY: all clean test

LIB_DIR:=lib
OBJ_DIR:=obj
BIN_DIR:=bin
SRC_DIR:=src
CPP_DIR:=cpp

CXX=g++
CXXFLAGS=-O3 -std=c++11 -fopenmp -g
OBJ=cpp/vg.pb.o xg.o # main.o not included for easier libxg.a creation
LD_INCLUDES=-I./ -Icpp -Istream/src -I$(SRC_DIR)
LD_LIBS=-lprotobuf -lsdsl -lz -ldivsufsort -ldivsufsort64 -lgomp -lm -lpthread
STREAM=stream
EXE:=xg


#Some little adjustments to build on OSX
#(tested with gcc4.9 installed from MacPorts)
SYS=$(shell uname -s)
all: $(BIN_DIR)/$(EXE) $(LIB_DIR)/libxg.a

doc: README.md
README.md: README.base.md
	pandoc -o README.html -s README.base.md
	pandoc -o DESIGN.html -s DESIGN.md
	cat README.base.md >README.md
	cat DESIGN.html | tail -7| perl -p -e 's/<p>/\n/g' | sed 's%</p>%%g' | head -10 >>README.md



$(CPP_DIR)/vg.pb.cc: $(CPP_DIR)/vg.pb.h
$(CPP_DIR)/vg.pb.h: $(SRC_DIR)/vg.proto
	mkdir -p cpp
	protoc $(SRC_DIR)/vg.proto --proto_path=$(SRC_DIR) --cpp_out=cpp

$(OBJ_DIR)/vg.pb.o: $(CPP_DIR)/vg.pb.h $(CPP_DIR)/vg.pb.cc
	$(CXX) $(CXXFLAGS) -c -o $(CPP_DIR)/vg.pb.o $(CPP_DIR)/vg.pb.cc $(LD_INCLUDES) $(LD_LIBS)

$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp $(CPP_DIR)/vg.pb.h $(SRC_DIR)/xg.hpp 
	$(CXX) $(CXXFLAGS) $(LD_LIBS) -c -o $@ $(SRC_DIR)/main.cpp $(LD_INCLUDES)

$(OBJ_DIR)/xg.o: $(SRC_DIR)/xg.cpp $(SRC_DIR)/xg.hpp $(CPP_DIR)/vg.pb.h
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDES) $(LD_LIBS)

$(BIN_DIR)/$(EXE): $(OBJ_DIR)/main.o $(CPP_DIR)/vg.pb.o $(OBJ_DIR)/xg.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_INCLUDES) $(LD_LIBS) $(STATICFLAGS)

$(LIB_DIR)/libxg.a: $(CPP_DIR)/vg.pb.o $(OBJ_DIR)/xg.o
	ar rs $@ $(OBJ_DIR)/xg.o $(CPP_DIR)/vg.pb.o

test:
	cd test && make

clean-xg:
	rm -f *.o xg cpp/* libxg.a

clean:
	rm -rf cpp
	rm -f $(EXECUTABLE)
	rm -f *.o
	rm -f libxg.a
