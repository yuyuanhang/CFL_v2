cmake_minimum_required(VERSION 2.8)
project(SubgraphMatching)

include_directories(
        ${PROJECT_SOURCE_DIR}/
)

set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -g -Wall -march=native -pthread")

add_executable(CFL main.cc)

add_subdirectory(graph)
add_subdirectory(matching)
add_subdirectory(utility)

target_link_libraries(CFL PUBLIC graph)
target_link_libraries(CFL PUBLIC matching)
target_link_libraries(CFL PUBLIC utility)