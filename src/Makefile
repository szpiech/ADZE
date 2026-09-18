# adze build
#
# Plain make; see CMakeLists.txt for the CMake build, which detects OpenMP and
# provides an install target.
#
#   make                      serial build
#   make OPENMP=1             GCC or Linux clang
#   make OPENMP=1 OMPFLAGS="-Xpreprocessor -fopenmp" OMPLIBS=-lomp
#                             Apple clang (needs libomp installed)
#   make ZLIB=1               compressed (.gz/.bgz) input, needs zlib
#   make clean
#
# Both are optional; nothing else is needed. 1.0 required GSL for one
# combination iterator, which adze now generates itself.

CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall
LDLIBS   ?= -lm
PREFIX   ?= /usr/local

OMPFLAGS ?= -fopenmp
OMPLIBS  ?= -fopenmp

ifdef OPENMP
CXXFLAGS += $(OMPFLAGS)
LDLIBS   += $(OMPLIBS)
endif

# zlib gives .vcf.gz and .stru.gz input; without it those files are refused
# with a message rather than misparsed.
ifdef ZLIB
CXXFLAGS += -DADZE_HAVE_ZLIB
LDLIBS   += -lz
endif

# Hooks for flags that belong to a particular build rather than to the
# program: a static link for a binary that has to run on a machine without the
# toolchain, or several -arch flags for a universal macOS binary. They exist
# because setting CXXFLAGS or LDFLAGS on the command line replaces the value
# in this file, and make then ignores the += above it -- so ZLIB=1 would
# silently stop defining ADZE_HAVE_ZLIB. Appending a command-line variable is
# safe, which is what these are.
#
#   make ZLIB=1 EXTRA_LDFLAGS=-static
#   make ZLIB=1 EXTRA_CXXFLAGS="-arch x86_64 -arch arm64" \
#               EXTRA_LDFLAGS="-arch x86_64 -arch arm64"
CXXFLAGS += $(EXTRA_CXXFLAGS)
LDFLAGS  += $(EXTRA_LDFLAGS)

# Set I_PATH/L_PATH if headers or libraries live outside the default paths.
ifdef I_PATH
CXXFLAGS += -I$(I_PATH)
endif
ifdef L_PATH
LDFLAGS  += -L$(L_PATH)
endif

# On Windows (MSYS2/MinGW) the compiler appends .exe to an output name that
# carries no suffix, so make never sees the target it just built and relinks
# on every invocation. Naming it explicitly avoids that:
#   make TARGET=adze.exe
TARGET  = adze
SOURCES = ADZE_main.cpp ADZE_main_tools.cpp ADZE_pop.cpp ADZE_pfile.cpp \
          ADZE_stats.cpp ADZE_pbar.cpp ADZE_binom.cpp
OBJECTS = $(SOURCES:.cpp=.o)
HEADERS = ADZE_main_tools.h ADZE_pop.h ADZE_pfile.h ADZE_stats.h \
          ADZE_pbar.h ADZE_binom.h

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

# Every object depends on every header. The headers are small and few, and 1.0's
# Makefile listed no header dependencies at all -- so editing a header left
# stale objects linked against a different struct layout.
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

install: $(TARGET)
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp $(TARGET) $(DESTDIR)$(PREFIX)/bin/

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all install clean
