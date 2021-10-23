cmake_minimum_required(VERSION 3.20)
project(gravity_simulator)

set(CMAKE_CXX_STANDARD 14)
add_compile_options(-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors)

add_executable(gravity_simulator sim_soa.cpp sim_soa.hpp)
