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
target_include_directories (Test PRIVATE
    ${PROJECT_SOURCE_DIR}/include

# external library
    ${EIGEN3_INCLUDE_DIR}
    ${json_SOURCE_DIR}/include
)

target_link_libraries(Test
    gtest_main

# external library
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
)

# CTestへ単体テストを登録
include(GoogleTest)
gtest_add_tests(TARGET Test)