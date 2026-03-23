cmake_minimum_required (VERSION 3.8)

# project ("steering_functions")
find_package(Eigen3 REQUIRED)

file(GLOB SRC_HC_CC_SS src/hc_cc_state_space/*.cpp)
set(SRC_FILES
    src/dubins_state_space/dubins_state_space.cpp
    src/reeds_shepp_state_space/reeds_shepp_state_space.cpp
    src/utilities/utilities.cpp
    src/filter/ekf.cpp
    ${SRC_HC_CC_SS})

add_library(steering_functions ${SRC_FILES})
target_include_directories(steering_functions PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}/include ${EIGEN3_INCLUDE_DIRS})

target_link_libraries(steering_functions PUBLIC 
  common_headers)

# add_executable (steering_functions_node "steering_functions.cpp" "steering_functions.h" ${SRC_FILES})
