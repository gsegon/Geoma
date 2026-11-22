add_executable(test_algorithms
        test_algorithms_normalize.cpp
        test_algorithms_canonicalize.cpp
        test_algorithms_verify_equal.cpp
        test_algorithms_find_transform.cpp

)

target_link_libraries(test_algorithms PRIVATE GTest::gtest_main geoma)