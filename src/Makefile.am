bin_PROGRAMS = pisa

# Includes
pisa_CPPFLAGS = \
	-I$(top_srcdir)/deps/SeqLib \
	-I$(top_srcdir)/deps/SeqLib/htslib \
	-I$(top_srcdir)/deps/CTPL \
	-I$(top_srcdir)/src -Wno-sign-compare

# Static libs
pisa_LDADD = \
	$(top_builddir)/deps/SeqLib/src/libseqlib.a \
	$(top_builddir)/deps/SeqLib/bwa/libbwa.a \
	$(top_builddir)/deps/SeqLib/htslib/libhts.a \
	$(top_builddir)/deps/SeqLib/fermi-lite/libfml.a

# Source files
pisa_SOURCES = \
	main.cpp \
	$(top_builddir)/src/readCorrection.cpp \
	$(top_builddir)/src/readFilter.cpp \
	$(top_builddir)/src/ThreadedFermiAssembler.cpp \
    $(top_builddir)/src/pisaExtract.cpp \
    $(top_builddir)/src/pisaAssemble.cpp \
    $(top_builddir)/src/pisaTest.cpp \
    $(top_builddir)/src/utility.cpp
