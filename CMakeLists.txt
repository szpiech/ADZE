# adze -- allelic richness and private allelic richness by rarefaction
#
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build -j
#   ctest --test-dir build            # differential test vs. a 1.0 reference
#   cmake --install build --prefix ~/.local
#
# OpenMP is used when the toolchain provides it (on Apple clang that means
# libomp, e.g. `brew install libomp`), and the program builds and runs fine
# without it. There are no other dependencies: ADZE 1.0 required GSL for a
# single combination iterator, which adze now generates itself.

cmake_minimum_required(VERSION 3.16)
project(adze VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# So an installed binary still finds a libomp that lives outside the default
# search paths (a Homebrew or conda prefix, say).
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(adze
  src/ADZE_main.cpp
  src/ADZE_main_tools.cpp
  src/ADZE_pop.cpp
  src/ADZE_pfile.cpp
  src/ADZE_stats.cpp
  src/ADZE_pbar.cpp
  src/ADZE_binom.cpp)

target_include_directories(adze PRIVATE src)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(adze PRIVATE -Wall)
endif()

option(ADZE_USE_OPENMP "Use OpenMP if available" ON)

if(ADZE_USE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(adze PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "adze: OpenMP enabled (--threads N)")
  else()
    message(STATUS "adze: no OpenMP found; building serial (--threads ignored)")
  endif()
endif()

install(TARGETS adze RUNTIME DESTINATION bin)

# The regression suite is differential: it needs a binary built from the 1.0
# sources to compare against. Point ADZE_REFERENCE at one to enable the test.
include(CTest)
set(ADZE_REFERENCE "" CACHE FILEPATH "Path to an ADZE 1.0 binary for regression tests")

if(BUILD_TESTING AND ADZE_REFERENCE)
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_FOUND)
    add_test(NAME regression
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/regress.py
                     --candidate $<TARGET_FILE:adze>
                     --reference ${ADZE_REFERENCE}
                     --workdir ${CMAKE_CURRENT_BINARY_DIR}/regress)
  endif()
endif()
