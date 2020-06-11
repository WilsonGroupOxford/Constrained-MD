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

include(GoogleTest)
# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# Now simply link against gtest or gtest_main as needed. Eg
list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/src/main.cpp")
file(GLOB TESTSOURCES "test/*.cpp")
add_executable(C-MD-Test  ${TESTSOURCES} ${SOURCES})
target_compile_options(C-MD-Test PRIVATE -Wall -Wextra -Wpedantic  -static -static-libstdc++)
target_include_directories(C-MD-Test PUBLIC "${PROJECT_SOURCE_DIR}/include")
set_target_properties(C-MD-Test  PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)
target_link_libraries(C-MD-Test gtest_main)
target_link_libraries(C-MD-Test Eigen3::Eigen)
gtest_add_tests(TARGET C-MD-Test)
enable_testing()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/modules)
include(CodeCoverage)
append_coverage_compiler_flags()
setup_target_for_coverage_gcovr_html(
     NAME ctest_coverage                    # New target name
     EXECUTABLE ctest -j ${PROCESSOR_COUNT} # Executable in PROJECT_BINARY_DIR
     DEPENDENCIES C-MD-Test      # Dependencies to build first
     EXCLUDE ./test/*.cpp ./build/*
)

