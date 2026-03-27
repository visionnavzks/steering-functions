find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)

find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -m pybind11 --cmakedir
        RESULT_VARIABLE pybind11_cmake_dir_result
        OUTPUT_VARIABLE pybind11_cmake_dir
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(pybind11_cmake_dir_result EQUAL 0 AND EXISTS "${pybind11_cmake_dir}/pybind11Config.cmake")
        list(PREPEND CMAKE_PREFIX_PATH "${pybind11_cmake_dir}")
        find_package(pybind11 CONFIG REQUIRED)
    else()
        message(FATAL_ERROR
            "pybind11 CMake package not found. Install pybind11 for ${Python3_EXECUTABLE} or set pybind11_DIR.")
    endif()
endif()

pybind11_add_module(steering_functions_python MODULE module.cpp)
target_link_libraries(steering_functions_python PRIVATE steering_path_lib)
set_target_properties(steering_functions_python PROPERTIES OUTPUT_NAME _core)

add_custom_command(
    TARGET steering_functions_python
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:steering_functions_python>
            "${CMAKE_SOURCE_DIR}/python/steering_functions/$<TARGET_FILE_NAME:steering_functions_python>"
    COMMENT "Copying Python extension into python/steering_functions"
)

install(TARGETS steering_functions_python
        LIBRARY DESTINATION steering_functions)