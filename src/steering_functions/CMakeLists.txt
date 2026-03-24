cmake_minimum_required (VERSION 3.8)

# project ("steering_functions")

option(STEERING_FUNCTIONS_ENABLE_EXTENDED_TESTS
  "Run the upstream steering_functions integration test executable in ctest"
  OFF)

file(GLOB SRC_HC_CC_SS src/hc_cc_state_space/*.cpp)
set(SRC_FILES
    src/base_state_space.cpp
    src/dubins_state_space/dubins_state_space.cpp
    src/reeds_shepp_state_space/reeds_shepp_state_space.cpp
    src/utilities/utilities.cpp
    ${SRC_HC_CC_SS})

add_library(steering_functions ${SRC_FILES})
target_compile_features(steering_functions PUBLIC cxx_std_17)
target_include_directories(steering_functions PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(steering_functions PUBLIC 
  steering_state)

if(TARGET common_headers)
  target_link_libraries(steering_functions PUBLIC common_headers)
endif()

if(BUILD_TESTING)
  find_package(GTest REQUIRED)
  include(GoogleTest)

  add_executable(fresnel_test test/fresnel_test.cpp)
  target_link_libraries(fresnel_test PRIVATE steering_functions GTest::gtest_main)
  gtest_discover_tests(fresnel_test)

  add_executable(hc_cc_circle_test test/hc_cc_circle_test.cpp)
  target_link_libraries(hc_cc_circle_test PRIVATE steering_functions GTest::gtest_main)
  gtest_discover_tests(hc_cc_circle_test)

  add_executable(steering_functions_test test/steering_functions_test.cpp)
  target_link_libraries(steering_functions_test PRIVATE steering_functions GTest::gtest_main)

  if(STEERING_FUNCTIONS_ENABLE_EXTENDED_TESTS)
    gtest_discover_tests(steering_functions_test)
  endif()
endif()

# add_executable (steering_functions_node "steering_functions.cpp" "steering_functions.h" ${SRC_FILES})
