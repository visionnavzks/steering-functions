cmake_minimum_required(VERSION 3.10)
project(steering_refactor LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(steering_refactor
  src/angle_utils.cpp
)

target_include_directories(steering_refactor
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

enable_testing()

add_executable(angle_utils_test
  tests/angle_utils_test.cpp
)

target_link_libraries(angle_utils_test
  PRIVATE
    steering_refactor
)

add_test(NAME angle_utils_test COMMAND angle_utils_test)
