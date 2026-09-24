cmake_minimum_required(VERSION 3.15)
project(bh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_library(gf2k STATIC
    src/GF2k.cpp
    src/Poly.cpp
)
target_include_directories(gf2k PUBLIC include)   # <-- was src, should be include

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(gf2k PRIVATE -march=native -O3)
elseif(MSVC)
    target_compile_options(gf2k PRIVATE /arch:AVX2 /O2)
endif()

enable_testing()
add_subdirectory(tests)