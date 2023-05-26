cmake_minimum_required(VERSION 3.25)

set(LIBRARY_NAME GaussQuadrature)
add_library(${LIBRARY_NAME} STATIC GaussQuadrature.cpp)

# includeパスに追加
target_include_directories(${LIBRARY_NAME} PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${json_SOURCE_DIR}/include
)

# ライブラリにリンク
target_link_libraries(${LIBRARY_NAME} PRIVATE 
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
)

# コンパイラ関連
target_compile_options(${LIBRARY_NAME} PUBLIC -O3 -march=native)