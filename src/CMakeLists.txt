add_library(geoma INTERFACE)
target_sources(geoma INTERFACE
        include/geoma/geoma.h
        include/geoma/algorithms.h
        include/geoma/config.h
        include/geoma/primitives.h
)

target_include_directories(geoma INTERFACE include)

if (BUILD_TESTING)
    include(CTest)

    include(FetchContent)
    FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.18.0.zip
    )

    # Required for MSVC on Windows
    if (MSVC)
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    endif ()

    FetchContent_MakeAvailable(googletest)


    add_subdirectory(test)
endif ()