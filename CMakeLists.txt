cmake_minimum_required(VERSION 3.22.1)
project(libtsp)

enable_testing()

set(CMAKE_CXX_STANDARD 20)

include(cmake/testing.cmake)

add_library(tsp SHARED src/tsp.cpp src/held_karp_tsp.cpp src/heuristic_tsp.cpp)
target_include_directories(tsp PUBLIC include)

add_subdirectory(third_party/googletest)

add_subdirectory(tests)
add_subdirectory(benchmarks)