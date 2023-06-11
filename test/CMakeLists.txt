cmake_minimum_required(VERSION 3.25)

# GoogleTestをダウンロード
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.12.0.zip
)

# 親プロジェクトのコンパイラ・リンカ設定を上書きするのを防ぐ（Windowsのみ）
if(WIN32)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()

# googletestをターゲットとして使えるようにする
FetchContent_MakeAvailable(googletest)

enable_testing()

# 単体テストの定義
add_executable(Test test.cpp)

find_package(PkgConfig)

pkg_check_modules(EIGEN3  REQUIRED Eigen3)
pkg_check_modules(LAPACK  REQUIRED lapack)
pkg_check_modules(LAPACKE REQUIRED lapacke)
pkg_check_modules(BLAS REQUIRED blas)

target_include_directories(Test PRIVATE
    PUBLIC ${PROJECT_SOURCE_DIR}/include

    PUBLIC ${EIGEN3_INCLUDE_DIRS}
    PUBLIC ${LAPACK_INCLUDE_DIRS}
    PUBLIC ${LAPACKE_INCLUDE_DIRS}
    PUBLIC ${BLAS_INCLUDE_DIRS}
)

target_link_directories(Test PRIVATE
    PUBLIC ${LAPACK_LIBRARY_DIRS}
    PUBLIC ${LAPACKE_LIBRARY_DIRS}
    PUBLIC ${BLAS_LIBRARY_DIRS}
)

target_link_libraries(Test PRIVATE
    gtest_main
# external libs
    PUBLIC ${EIGEN3_LIBRARIES}
    PUBLIC ${LAPACK_LIBRARIES}
    PUBLIC ${LAPACKE_LIBRARIES}
    PUBLIC ${BLAS_LIBRARIES}
)

# CTestへ単体テストを登録
include(GoogleTest)
gtest_add_tests(TARGET Test)