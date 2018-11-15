AM_CXXFLAGS = @OPENMP_CXXFLAGS@ @PTHREAD_CFLAGS@ @PROTOBUF_CFLAGS@ @SDSL_CFLAGS@ @SPARSEHASH_CFLAGS@ @HTSLIB_CFLAGS@ @JANSSON_CFLAGS@

lib_LIBRARIES = libxg.a
libxg_a_SOURCES = xg.cpp xg.hpp hash_map.hpp stream.hpp stream.cpp alignment.cpp alignment.hpp edit.cpp edit.hpp json2pb.cpp json2pb.h types.hpp nodetraversal.hpp region.hpp region.cpp path.hpp path.cpp graph.hpp graph.cpp sha1.hpp sha1.cpp position.cpp position.hpp handle.hpp handle.cpp utility.hpp utility.cpp blocked_gzip_input_stream.cpp blocked_gzip_input_stream.hpp blocked_gzip_output_stream.cpp blocked_gzip_output_stream.hpp hfile_cppstream.hpp hfile_cppstream.cpp hfile_internal.h textutils_internal.h vg.pb.cc vg.pb.h

bin_PROGRAMS = xg
xg_SOURCES = main.cpp
xg_LDADD = libxg.a @PTHREAD_LIBS@ @PROTOBUF_LIBS@ @SDSL_LIBS@ @SPARSEHASH_LIBS@ @HTSLIB_LIBS@ @JANSSON_LIBS@
xg_LDFLAGS = @OPENMP_CXXFLAGS@

BUILT_SOURCES = hash_map.hpp stream.hpp edit.hpp json2pb.h types.hpp nodetraversal.hpp region.hpp graph.hpp blocked_gzip_input_stream.hpp blocked_gzip_output_stream.hpp hfile_cppstream.hpp sha1.hpp sparsepp hfile_internal.h textutils_internal.h vg.pb.h

hash_map.hpp:
	git submodule init && git submodule update
	cp ../vg/src/hash_map.hpp .

stream.cpp: stream.hpp
stream.hpp:
	cp ../vg/src/stream.hpp .
	cp ../vg/src/stream.cpp .

edit.cpp: edit.hpp
edit.hpp:
	cp ../vg/src/edit.hpp .
	cp ../vg/src/edit.cpp .

json2pb.cpp: json2pb.h
json2pb.h:
	cp ../vg/src/json2pb.h .
	cp ../vg/src/json2pb.cpp .

types.hpp:
	cp ../vg/src/types.hpp .

nodetraversal.hpp:
	cp ../vg/src/nodetraversal.hpp .

region.cpp: region.hpp
region.hpp:
	cp ../vg/src/region.hpp .
	cp ../vg/src/region.cpp .

graph.cpp: graph.hpp
graph.hpp:
	cp ../vg/src/graph.hpp .
	cp ../vg/src/graph.cpp .

blocked_gzip_input_stream.cpp: blocked_gzip_input_stream.hpp
blocked_gzip_input_stream.hpp:
	cp ../vg/src/blocked_gzip_input_stream.hpp .
	cp ../vg/src/blocked_gzip_input_stream.cpp .

blocked_gzip_output_stream.cpp: blocked_gzip_output_stream.hpp
blocked_gzip_output_stream.hpp:
	cp ../vg/src/blocked_gzip_output_stream.hpp .
	cp ../vg/src/blocked_gzip_output_stream.cpp .

hfile_cppstream.cpp: hfile_cppstream.hpp
hfile_cppstream.hpp:
	cp ../vg/src/hfile_cppstream.hpp .
	cp ../vg/src/hfile_cppstream.cpp .

sha1.cpp: sha1.hpp
sha1.hpp:
	cp ../sha1/sha1.hpp .
	cp ../sha1/sha1.cpp .

sparsepp:
	cp -r ../sparsepp/sparsepp .

hfile_internal.h:
	cp ../htslib/hfile_internal.h .

textutils_internal.h:
	cp ../htslib/textutils_internal.h .

vg.proto:
	cp ../vg/src/vg.proto .

vg.pb.cc: vg.pb.h
vg.pb.h: vg.proto
	protoc ./vg.proto --proto_path=./ --cpp_out=./
