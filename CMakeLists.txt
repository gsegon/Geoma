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

# executable (remove)
add_executable(untitled1 main.cpp)
target_link_libraries(untitled1 PUBLIC geoma)
