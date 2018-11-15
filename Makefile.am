AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src
dist_doc_DATA = README.md

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = xg.pc

pkginclude_HEADERS = src/xg.hpp src/hash_map_set.hpp src/vg.pb.h

DISTCLEANFILES = src/stream.hpp src/vg.pb.h src/vg.pb.cc

doc: README.md
README.md: README.base.md
	pandoc -o README.html -s README.base.md
	pandoc -o DESIGN.html -s DESIGN.md
	cat README.base.md >README.md
	cat DESIGN.html | tail -7| perl -p -e 's/<p>/\n/g' | sed 's%</p>%%g' | head -10 >>README.md
