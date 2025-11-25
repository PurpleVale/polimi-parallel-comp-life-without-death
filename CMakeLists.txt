cmake_minimum_required(VERSION 3.28)
set(CMAKE_CXX_STANDARD 20)

project(Life_Without_Death)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP Correctly Initialized")
endif()

set(DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${DIR}/bin")

add_executable(Base src/base_implementation.cpp)
target_link_libraries(Base PUBLIC OpenMP::OpenMP_CXX)

add_executable(LwD src/our_implementation.cpp)
target_link_libraries(LwD PUBLIC OpenMP::OpenMP_CXX)