cmake_minimum_required(VERSION 3.16)

project(steering_functions_workspace LANGUAGES CXX)

include(CTest)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(steering_state INTERFACE)
target_compile_features(steering_state INTERFACE cxx_std_17)
target_include_directories(steering_state
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

add_subdirectory(src/steering_functions)
add_subdirectory(src/steering_path_lib)