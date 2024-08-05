cmake_minimum_required(VERSION 3.30)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# project name
project(dev)

# set c++-17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set export compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set compiler flags
message(STATUS ${CMAKE_CXX_COMPILER_ID})

# すべてのライブラリを静的にビルド
set(BUILD_SHARED_LIBS OFF)

# すべてのライブラリを静的にビルド
set(BUILD_SHARED_LIBS OFF)

# リリースモードの最適化
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Ot")

# Windowsの標準ライブラリの静的リンク
if(WIN32)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")

    # MSVCでOpenMPを有効にする
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
endif()

# download external libraries
include(ExternalProject)

# download eigen library
set(EIGEN_BUILD_DIR ${CMAKE_BINARY_DIR}/eigen)
set(EIGEN_INSTALL_DIR ${PROJECT_SOURCE_DIR}/include/eigen3)
set(EIGEN_INCLUDE_DIR ${EIGEN_INSTALL_DIR})

if(NOT OFFLINE)
    ExternalProject_Add(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
        PREFIX ${EIGEN_BUILD_DIR}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        UPDATE_COMMAND ""
        INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy_directory ${EIGEN_BUILD_DIR}/src/eigen/Eigen ${EIGEN_INCLUDE_DIR}/Eigen
        && ${CMAKE_COMMAND} -E copy_directory ${EIGEN_BUILD_DIR}/src/eigen/unsupported ${EIGEN_INCLUDE_DIR}/unsupported
        TEST_COMMAND ""
    )
endif()

include_directories(SYSTEM ${EIGEN_INCLUDE_DIR})

# main executable file list
# search .cpp files
file(GLOB SOURCES
    ./src/*.cpp)

# add all .cpp files to target
add_executable(exe ${SOURCES})

# process before build
if(NOT OFFLINE)
    add_dependencies(exe eigen)
endif()