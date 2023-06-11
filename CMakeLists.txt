cmake_minimum_required(VERSION 3.25)

set(PROJECT_NAME main)
project (${PROJECT_NAME})

add_executable(${PROJECT_NAME} main.cpp)

find_package(PkgConfig)

pkg_check_modules(EIGEN3  REQUIRED Eigen3)
pkg_check_modules(LAPACK  REQUIRED lapack)
pkg_check_modules(LAPACKE REQUIRED lapacke)
pkg_check_modules(BLAS REQUIRED blas)

target_include_directories(${PROJECT_NAME} 
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    PUBLIC ${EIGEN3_INCLUDE_DIRS}
    PUBLIC ${LAPACK_INCLUDE_DIRS}
    PUBLIC ${LAPACKE_INCLUDE_DIRS}
    PUBLIC ${BLAS_INCLUDE_DIRS}
)

target_link_directories(${PROJECT_NAME} 
    PUBLIC ${LAPACK_LIBRARY_DIRS}
    PUBLIC ${LAPACKE_LIBRARY_DIRS}
    PUBLIC ${BLAS_LIBRARY_DIRS}
)

add_subdirectory(test)

target_link_libraries(${PROJECT_NAME} PRIVATE
# external libs
    PUBLIC ${EIGEN3_LIBRARIES}
    PUBLIC ${LAPACK_LIBRARIES}
    PUBLIC ${LAPACKE_LIBRARIES}
    PUBLIC ${BLAS_LIBRARIES}
)

# コンパイラ設定
set(CMAKE_C_COMPILER /usr/local/bin/gcc-13)
set(CMAKE_CXX_COMPILER /usr/local/bin/g++-13)

# コンパイラ関連
target_compile_options(${PROJECT_NAME} PUBLIC -O3 -march=native)