cmake_minimum_required(VERSION 3.20)
project(gravity_simulator)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(gravity_simulator
        sim-v3
        sim_soa.cpp
        sim_soa.hpp)
