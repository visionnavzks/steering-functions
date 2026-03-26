cmake_minimum_required(VERSION 3.12)
project(steering_path_lib)

# Create the library
add_library(steering_path_lib)
target_compile_features(steering_path_lib PUBLIC cxx_std_17)

# Specify the include directories for the library
target_include_directories(steering_path_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Add source files
target_sources(steering_path_lib
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/steering_path.cpp
)

# Link against required dependencies
target_link_libraries(steering_path_lib
    PUBLIC
        steering_functions
)
