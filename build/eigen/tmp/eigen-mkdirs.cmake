# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "C:/Users/Shimolab30/OneDrive/c++/minimalFEM-VTK/build/eigen/src/eigen")
  file(MAKE_DIRECTORY "C:/Users/Shimolab30/OneDrive/c++/minimalFEM-VTK/build/eigen/src/eigen")
endif()
file(MAKE_DIRECTORY
  "C:/Users/Shimolab30/OneDrive/c++/minimalFEM-VTK/build/eigen/src/eigen-build"
  "C:/Users/Shimolab30/OneDrive/c++/minimalFEM-VTK/build/eigen"
  "C:/Users/Shimolab30/OneDrive/c++/minimalFEM-VTK/build/eigen/tmp"
  "C:/Users/Shimolab30/OneDrive/c++/minimalFEM-VTK/build/eigen/src/eigen-stamp"
  "C:/Users/Shimolab30/OneDrive/c++/minimalFEM-VTK/build/eigen/src"
  "C:/Users/Shimolab30/OneDrive/c++/minimalFEM-VTK/build/eigen/src/eigen-stamp"
)

set(configSubDirs Debug;Release;MinSizeRel;RelWithDebInfo)
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "C:/Users/Shimolab30/OneDrive/c++/minimalFEM-VTK/build/eigen/src/eigen-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "C:/Users/Shimolab30/OneDrive/c++/minimalFEM-VTK/build/eigen/src/eigen-stamp${cfgdir}") # cfgdir has leading slash
endif()
