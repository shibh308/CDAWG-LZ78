cmake_minimum_required(VERSION 3.22)
project(CDAWG_LZ78)

set(CMAKE_CXX_STANDARD 20)

set(SDSL_INCLUDE_DIR "./sdsl-lite/include")
set(SDSL_LIBRARY_DIR "./sdsl-lite/lib")

include_directories(${SDSL_INCLUDE_DIR})
link_directories(${SDSL_LIBRARY_DIR})

add_executable(CDAWG_LZ78 main.cpp
        includes/map.hpp
        includes/vector.hpp
        includes/dawg.hpp
        includes/suffix_array.hpp
        includes/lz78.hpp
        includes/markedancestor.hpp
        includes/biasedsearchtree.hpp
        includes/sdsl_suffixtree.hpp
)

target_link_libraries(CDAWG_LZ78 divsufsort divsufsort64 sdsl)