cmake_minimum_required(VERSION 3.22.1)
project(libtsp)

enable_testing()

set(CMAKE_CXX_STANDARD 20)

include(cmake/testing.cmake)

option(LIBTSP_BUILD_STATIC "Build libtsp as a static library" ON)

set(LIBTSP_SOURCES src/tsp.cpp src/held_karp_tsp.cpp src/heuristic_tsp.cpp)

if (LIBTSP_BUILD_STATIC)
    add_library(tsp STATIC ${LIBTSP_SOURCES})
else()
    add_library(tsp SHARED ${LIBTSP_SOURCES})
endif ()

target_include_directories(tsp PUBLIC include)

add_subdirectory(third_party/googletest)

add_subdirectory(tests)
add_subdirectory(benchmarks)