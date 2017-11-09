project(sLLGS C CXX Fortran)
cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 14)



# gfortran
#set (CMAKE_Fortran_FLAGS_RELEASE "-ffree-line-length-none -funroll-all-loops -fno-f2c -O3")
#set (CMAKE_Fortran_FLAGS_DEBUG   "-ffree-line-length-none -fno-f2c -O0 -g")

aux_source_directory(. SRC_LIST)

set(FORTRAN_SRC demag3D3_modificado.f90 dipolar3D3_modificado.f90)
add_executable(${PROJECT_NAME} ${SRC_LIST} ${FORTRAN_SRC})

## To use MatplotLib
find_package(PythonLibs 2.7)
target_include_directories(${PROJECT_NAME} PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES})