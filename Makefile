.PHONY: all clean test pre

LIB_DIR:=lib
INC_DIR:=include
OBJ_DIR:=obj
BIN_DIR:=bin
SRC_DIR:=src
CPP_DIR:=cpp

CXX=g++
CXXFLAGS=-O3 -std=c++11 -fopenmp -ggdb
OBJ=cpp/vg.pb.o xg.o # main.o not included for easier libxg.a creation
LD_INCLUDES=-I./ -Icpp -Istream/src -IDYNAMIC/include -IDYNAMIC/include/internal -IDYNAMIC/include/algorithms -I$(SRC_DIR)
LD_LIBS=-lprotobuf -lsdsl -lz -ldivsufsort -ldivsufsort64 -lgomp -lm -lpthread
STREAM=stream
EXE:=xg
CWD:=$(shell pwd)

#Some little adjustments to build on OSX
#(tested with gcc4.9 installed from MacPorts)
SYS=$(shell uname -s)
ifeq (${SYS},Darwin)
        CMAKE_BIN= # linux binary won't work on os x
        CMAKE_SETPATH= # mac ports cmake seems to work just fine so leave blank
        STATICFLAGS= # -static doesn't work on OSX unless libgcc compiled as static.
else
        CMAKE_BIN=cmake-3.3.0-rc2-Linux-x86_64/bin/cmake
        CMAKE_SETPATH=PATH=./cmake-3.3.0-rc2-Linux-x86_64/bin/:${PATH}
        STATICFLAGS=-static -static-libstdc++ -static-libgcc -Wl,-Bstatic
endif

all: $(BIN_DIR)/$(EXE) $(LIB_DIR)/libxg.a

doc: README.md
README.md: README.base.md
	pandoc -o README.html -s README.base.md
	pandoc -o DESIGN.html -s DESIGN.md
	cat README.base.md >README.md
	cat DESIGN.html | tail -7| perl -p -e 's/<p>/\n/g' | sed 's%</p>%%g' | head -10 >>README.md

pre:
	if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	if [ ! -d $(INC_DIR) ]; then mkdir -p $(INC_DIR); fi
	if [ ! -d $(CPP_DIR) ]; then mkdir -p $(CPP_DIR); fi

$(CMAKE_BIN):
	wget --no-check-certificate http://www.cmake.org/files/v3.3/cmake-3.3.0-rc2-Linux-x86_64.tar.gz
	tar xzvf cmake-3.3.0-rc2-Linux-x86_64.tar.gz

protobuf/:
	git clone https://github.com/google/protobuf.git

sdsl-lite/: $(CMAKE_BIN)
	git clone https://github.com/simongog/sdsl-lite.git

DYNAMIC/:
	git clone --recursive https://github.com/vgteam/DYNAMIC.git
	cd DYNAMIC && git checkout cfd3ae7

get-deps: $(CMAKE_BIN) protobuf/ sdsl-lite/ DYNAMIC/
	cd protobuf && git checkout dfae9e3 && ./autogen.sh || ./autogen.sh && ./configure --prefix="$(CWD)" && make -j 8 && make install && export PATH=$(CWD)/bin:$$PATH
	PATH=`pwd`/cmake-3.3.0-rc2-Linux-x86_64/bin/:$$PATH && cd sdsl-lite && git checkout 25b20b0 && ./install.sh $(CWD)

$(CPP_DIR)/vg.pb.cc: $(CPP_DIR)/vg.pb.h
$(CPP_DIR)/vg.pb.h: $(SRC_DIR)/vg.proto pre
	mkdir -p cpp
	protoc $(SRC_DIR)/vg.proto --proto_path=$(SRC_DIR) --cpp_out=cpp

$(OBJ_DIR)/vg.pb.o: $(CPP_DIR)/vg.pb.h $(CPP_DIR)/vg.pb.cc pre
	$(CXX) $(CXXFLAGS) -c -o $(CPP_DIR)/vg.pb.o $(CPP_DIR)/vg.pb.cc $(LD_INCLUDES) $(LD_LIBS)

$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp $(CPP_DIR)/vg.pb.h $(SRC_DIR)/xg.hpp pre
	$(CXX) $(CXXFLAGS) $(LD_LIBS) -c -o $@ $(SRC_DIR)/main.cpp $(LD_INCLUDES)

$(OBJ_DIR)/xg.o: $(SRC_DIR)/xg.cpp $(SRC_DIR)/xg.hpp $(CPP_DIR)/vg.pb.h pre
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDES) $(LD_LIBS)

$(BIN_DIR)/$(EXE): $(OBJ_DIR)/main.o $(CPP_DIR)/vg.pb.o $(OBJ_DIR)/xg.o $(INC_DIR)/stream.h pre 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ_DIR)/main.o $(CPP_DIR)/vg.pb.o $(OBJ_DIR)/xg.o $(LD_INCLUDES) $(LD_LIBS) $(STATICFLAGS)

$(LIB_DIR)/libxg.a: $(CPP_DIR)/vg.pb.o $(OBJ_DIR)/xg.o $(INC_DIR)/stream.h pre
	ar rs $@ $(OBJ_DIR)/xg.o $(CPP_DIR)/vg.pb.o

$(INC_DIR)/stream.h: pre 
	cd stream && $(MAKE) && cp include/* ../include/

test:
	cd test && make

clean-xg:
	rm -f *.o xg cpp/* libxg.a

clean:
	rm -rf cpp
	rm -f $(EXECUTABLE)
	rm -f *.o
	rm -f libxg.a
	rm -rf $(INC_DIR)
	rm -rf $(LIB_DIR)
	rm -rf $(BIN_DIR)
	rm -rf $(OBJ_DIR)
