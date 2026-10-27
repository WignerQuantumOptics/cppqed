# Copyright András Vukics 2020–2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDutils")

# ### Library detection
find_package(Boost ${BOOST_MINIMAL} REQUIRED COMPONENTS iostreams json serialization)

find_package(Eigen3 3.4 REQUIRED)

# ### Source dependencies, pinned to full commit SHAs
#
# Offline builds (cluster compute nodes): -DFETCHCONTENT_FULLY_DISCONNECTED=ON together with
# -DFETCHCONTENT_SOURCE_DIR_<NAME>=<local copy> for each of the names below (upper case).
# Co-development against a sibling library checkout:
# -DFETCHCONTENT_SOURCE_DIR_CPPQED-MULTIDIAGONAL=../cppqed-multidiagonal
include(FetchContent)

# An installed library is used when present; 0.4 is the first version with twoITimesImagPartOf
# (headers are in include/multidiagonal/ since 0.3).
FetchContent_Declare(cppqed-multidiagonal
  GIT_REPOSITORY    https://github.com/WignerQuantumOptics/cppqed-multidiagonal.git
  GIT_TAG           24d651cddade542fa2cb193f3c984463385d4017 # v0.4.0, 2026-10-27
  FIND_PACKAGE_ARGS 0.4)

# header-only; popl's own CMakeLists (examples, tests) is skipped by pointing SOURCE_SUBDIR at include/
FetchContent_Declare(popl
  GIT_REPOSITORY https://github.com/badaix/popl.git
  GIT_TAG        bda5f43099d67419089a44c9e54474e4998a9a26  # master, 2021-08-02
  SOURCE_SUBDIR  include)
FetchContent_Declare(pcg-cpp
  GIT_REPOSITORY https://github.com/imneme/pcg-cpp.git
  GIT_TAG        428802d1a5634f96bcd0705fab379ff0113bcf13) # master, 2022-04-08
FetchContent_Declare(xoshiro-cpp
  GIT_REPOSITORY https://github.com/Reputeless/Xoshiro-cpp.git
  GIT_TAG        19bcbb2ce0ed158233187f524fd0964c105a65b3) # master, 2021-08-05

FetchContent_MakeAvailable(cppqed-multidiagonal popl pcg-cpp xoshiro-cpp)

add_library(ThirdParty INTERFACE)
target_include_directories(ThirdParty INTERFACE ${xoshiro-cpp_SOURCE_DIR} ${pcg-cpp_SOURCE_DIR}/include ${popl_SOURCE_DIR}/include)

add_library(CPPQEDutils SHARED
            $<TARGET_OBJECTS:cppqedGitSHA1>
            Pars.cc Trajectory.cc Version.cc
            )

add_library(CPPQED::CPPQEDutils ALIAS CPPQEDutils)

target_link_libraries(CPPQEDutils PUBLIC ThirdParty cppqed::multidiagonal Boost::json Boost::serialization Eigen3::Eigen)

install(TARGETS CPPQEDutils DESTINATION lib)