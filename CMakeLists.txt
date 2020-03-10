cmake_minimum_required(VERSION 3.5)

project(Constrained-MD LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


file(GLOB SOURCES "src/*.cpp")
add_executable(Constrained-MD ${SOURCES})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(Constrained-MD Eigen3::Eigen)

target_compile_options(Constrained-MD PRIVATE -Wall -Wextra -Wpedantic)
target_include_directories(Constrained-MD PUBLIC "${PROJECT_SOURCE_DIR}/include")
