add_library(geoma INTERFACE)
target_sources(geoma INTERFACE
        include/geoma/geoma.h
        include/geoma/algorithms.h
        include/geoma/config.h
        include/geoma/primitives.h
)

target_include_directories(geoma INTERFACE include)

if(BUILD_TESTING)
    include(CTest)

    FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.zip
    )

    # Required for MSVC on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(googletest)


    add_subdirectory(test)
endif ()