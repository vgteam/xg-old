AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src
dist_doc_DATA = README.md

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = xg.pc

pkginclude_HEADERS = src/xg.hpp src/hash_map.hpp src/stream.hpp src/alignment.hpp src/edit.hpp src/json2pb.h src/types.hpp src/nodetraversal.hpp src/region.hpp src/path.hpp src/graph.hpp src/sha1.hpp src/position.hpp src/utility.hpp src/handle.hpp src/blocked_gzip_input_stream.hpp src/blocked_gzip_output_stream.hpp src/hfile_cppstream.hpp src/hfile_internal.h src/textutils_internal.h src/vg.pb.h
sparseppincludedir = $(pkgincludedir)/sparsepp
sparseppinclude_HEADERS = src/sparsepp/spp_config.h  src/sparsepp/spp_dlalloc.h  src/sparsepp/spp.h  src/sparsepp/spp_memory.h  src/sparsepp/spp_smartptr.h  src/sparsepp/spp_stdint.h  src/sparsepp/spp_timer.h  src/sparsepp/spp_traits.h  src/sparsepp/spp_utils.h

DISTCLEANFILES = src/stream.hpp src/stream.cpp src/edit.hpp src/edit.cpp src/json2pb.h src/json2pb.cpp src/types.hpp src/nodetraversal.hpp src/region.hpp src/region.cpp src/graph.hpp src/graph.cpp src/blocked_gzip_input_stream.hpp src/blocked_gzip_input_stream.cpp src/blocked_gzip_output_stream.hpp src/blocked_gzip_output_stream.cpp src/hfile_cppstream.hpp src/hfile_cppstream.cpp src/sha1.hpp src/sha1.cpp src/vg.pb.h src/vg.pb.cc src/hfile_internal.h src/textutils_internal.h src/vg.proto

distclean-local:
	-rm -rf src/sparsepp

doc: README.md
README.md: README.base.md
	pandoc -o README.html -s README.base.md
	pandoc -o DESIGN.html -s DESIGN.md
	cat README.base.md >README.md
	cat DESIGN.html | tail -7| perl -p -e 's/<p>/\n/g' | sed 's%</p>%%g' | head -10 >>README.md
