cmake_minimum_required(VERSION 4.3)
project(Geoma VERSION 0.1.0 LANGUAGES CXX)

# C++ standard requirement
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Build options
option(BUILD_TESTING "Enable building tests" OFF)

# Add source
add_subdirectory(src)

# Add examples
add_subdirectory(examples)

