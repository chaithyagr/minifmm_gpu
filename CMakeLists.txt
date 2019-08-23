# CMakeLists.txt
cmake_minimum_required(VERSION 3.5)
project(minifmm)

find_package(PythonInterp 3.5 REQUIRED)
find_package(PythonLibs 3.5 REQUIRED)
add_subdirectory(pybind11)

include_directories(${Boost_INCLUDE_DIRS})
set(CUDA_SEPARABLE_COMPILATION ON)

find_package(CUDA REQUIRED)

# Gather all source files
FILE(GLOB sources src/*.c*)

include_directories(${CUDA_INCLUDE_DIRS})
include_directories(include/
        pybind11/include
)
include_directories(
    ${PYTHON_INCLUDE_DIRS}
    )
#set_source_files_properties( ${sources} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
#set_source_files_properties(main.cpp PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
set(CMAKE_CXX_FLAGS "-std=c++11 -fPIC")

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_52 -O3 --use_fast_math -D_FORCE_INLINES
        --compiler-options -fPIC --relocatable-device-code=true )

cuda_add_library(fmmlib SHARED ${sources} )

pybind11_add_module(minifmm main.cpp)
target_link_libraries(minifmm PRIVATE fmmlib)

#This is for debug of C++ code
add_executable(fmm run_cpp.cpp)
target_link_libraries(fmm fmmlib )
