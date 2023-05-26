cmake_minimum_required(VERSION 3.25)

# To use the googletest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Generate executable file
add_executable(Test ./test.cpp)

target_link_libraries(Test PRIVATE 
    GTest::gtest_main
    GaussQuadrature
)

include(GoogleTest)
gtest_discover_tests(Test)