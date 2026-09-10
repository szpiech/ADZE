# adze build
#
# Plain make; see CMakeLists.txt for the CMake build, which detects OpenMP and
# provides an install target.
#
#   make                      serial build
#   make OPENMP=1             GCC or Linux clang
#   make OPENMP=1 OMPFLAGS="-Xpreprocessor -fopenmp" OMPLIBS=-lomp
#                             Apple clang (needs libomp installed)
#   make clean
#
# There is no external dependency: 1.0 needed GSL for one combination iterator.

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

# Set I_PATH/L_PATH if headers or libraries live outside the default paths.
ifdef I_PATH
CXXFLAGS += -I$(I_PATH)
endif
ifdef L_PATH
LDFLAGS  += -L$(L_PATH)
endif

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
