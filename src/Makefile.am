AM_CXXFLAGS = @OPENMP_CXXFLAGS@ @PROTOBUF_CFLAGS@ @SDSL_CFLAGS@ @PTHREAD_CFLAGS@

lib_LIBRARIES = libxg.a
libxg_a_SOURCES = xg.cpp vg.pb.cc hash_map_set.hpp xg.hpp

bin_PROGRAMS = xg
xg_SOURCES = main.cpp
xg_LDADD = libxg.a @PROTOBUF_LIBS@ @PTHREAD_LIBS@ @SDSL_LIBS@
xg_LDFLAGS = @OPENMP_CXXFLAGS@

BUILT_SOURCES = stream.hpp vg.pb.h

stream.hpp:
	git submodule init && git submodule update
	cp ../stream/src/stream.hpp .

vg.pb.cc: vg.pb.h
vg.pb.h: vg.proto
	protoc ./vg.proto --proto_path=./ --cpp_out=./
