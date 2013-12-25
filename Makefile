# Makefile for SoftFM

# ----- Tweak these settings to configure for your system

# TODO : -D_FILE_OFFSET_BITS=64

CROSS         =
CFLAGS_OPT    = -O2 -ffast-math -ftree-vectorize
CFLAGS_DEBUG  = -g
CFLAGS_ARCH   =
CFLAGS_PATH   = -I/home/joris/test/rtl-sdr/inst/include
CFLAGS_EXTRA  =
LDFLAGS_PATH  = -L/home/joris/test/rtl-sdr/inst/lib
LDFLAGS_EXTRA =
LIBS_RTLSDR   = /home/joris/test/rtl-sdr/inst/lib/librtlsdr.a -lusb-1.0
LIBS_EXTRA    =

# ----- end tweakable settings


CXX = $(CROSS)g++
CXXFLAGS = -std=c++11 -Wall $(CFLAGS_OPT) $(CFLAGS_DEBUG) \
           $(CFLAGS_ARCH) $(CFLAGS_PATH) $(CFLAGS_EXTRA)
LDFLAGS = $(LDFLAGS_PATH) $(LDFLAGS_EXTRA)
LDLIBS  = $(LIBS_RTLSDR) $(LIBS_EXTRA)

OBJS	= RtlSdrSource.o Filter.o FmDecode.o AudioOutput.o main.o

softfm         : $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

RtlSdrSource.o : RtlSdrSource.cc RtlSdrSource.h SoftFM.h
Filter.o       : Filter.cc Filter.h SoftFM.h
FmDecode.o     : FmDecode.cc FmDecode.h SoftFM.h Filter.h
AudioOutput.o  : AudioOutput.cc AudioOutput.h SoftFM.h
main.o         : main.cc SoftFM.h RtlSdrSource.h Filter.h FmDecode.h AudioOutput.h

.PHONY: clean
clean:
	$(RM) softfm *.o

