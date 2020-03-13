cmake_minimum_required(VERSION 3.5)

project(Constrained-MD LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


file(GLOB SOURCES "src/*.cpp")
add_executable(Constrained-MD ${SOURCES})

target_include_directories(Constrained-MD PUBLIC "${PROJECT_SOURCE_DIR}/include")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(Constrained-MD Eigen3::Eigen)

find_package (yaml-cpp REQUIRED NO_MODULE)
target_link_libraries(Constrained-MD yaml-cpp)

target_compile_options(Constrained-MD PRIVATE -Wall -Wextra -Wpedantic)

