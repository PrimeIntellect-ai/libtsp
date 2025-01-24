add_executable(test_basic test_basic.cpp)
target_link_libraries(test_basic PRIVATE tsp gtest_main)

add_executable(test_complex test_complex.cpp)
target_link_libraries(test_complex PRIVATE tsp gtest_main)

add_subdirectory(test_data_gen)