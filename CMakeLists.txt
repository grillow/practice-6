project(safety_calculation)
cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 20)

add_executable(${PROJECT_NAME}
        src/main.cpp
        src/gps.hpp src/DataReader.hpp src/data.hpp)

find_package(HDF5 REQUIRED COMPONENTS HL)
target_link_libraries(${PROJECT_NAME} PRIVATE HDF5::HDF5)
