cmake_minimum_required(VERSION 3.22)
project(CDAWG_LZ78)

set(CMAKE_CXX_STANDARD 20)

add_executable(CDAWG_LZ78 main.cpp
        includes/map.hpp
        includes/vector.hpp
        includes/dawg.hpp
        includes/suffix_array.hpp
        includes/lz78.hpp
        includes/markedancestor.hpp
        includes/biasedsearchtree.hpp
)
