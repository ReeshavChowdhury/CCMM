cmake_minimum_required(VERSION 3.10)
project(AFCMT_Benchmarks)
set(CMAKE_CXX_STANDARD 17)

# Force the local install directory
set(OpenFHE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../install/lib/OpenFHE" CACHE PATH "OpenFHE config directory" FORCE)

find_package(OpenFHE REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(${OpenFHE_INCLUDE})
include_directories(${OpenFHE_INCLUDE}/core)
include_directories(${OpenFHE_INCLUDE}/pke)
include_directories(${OpenFHE_INCLUDE}/binfhe)
include_directories(${OpenFHE_INCLUDE}/third-party/include)
link_directories(${OpenFHE_LIBDIR})

add_executable(levelA_B_bench levelA_B_bench.cpp)
target_link_libraries(levelA_B_bench ${OpenFHE_LIBRARIES} OpenMP::OpenMP_CXX)

add_executable(rect_bench rect_bench.cpp)
target_link_libraries(rect_bench ${OpenFHE_LIBRARIES} OpenMP::OpenMP_CXX)

# Force RPATH to use the local install
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/../install/lib")
set_target_properties(levelA_B_bench PROPERTIES INSTALL_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/../install/lib")
set_target_properties(rect_bench PROPERTIES INSTALL_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/../install/lib")
