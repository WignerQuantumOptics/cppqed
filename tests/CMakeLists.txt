# Copyright András Vukics 2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring tests")

# Unit tests (Catch2 v3, as the library). An installed Catch2 3 is used when present.
include(FetchContent)
FetchContent_Declare(Catch2
  GIT_REPOSITORY    https://github.com/catchorg/Catch2.git
  GIT_TAG           fa43b77429ba76c462b1898d6cd2f2d7a9416b14 # v3.7.1
  FIND_PACKAGE_ARGS 3)
FetchContent_MakeAvailable(Catch2)

# testHeaders.cc goes in twice, as two translation units, to catch ODR violations in the headers (#51)
add_library(cppqedTestsHeadersCopy OBJECT unit/testHeaders.cc)
target_compile_definitions(cppqedTestsHeadersCopy PRIVATE ODR_COPY="2")
target_link_libraries(cppqedTestsHeadersCopy PRIVATE CPPQED::CPPQEDelements Catch2::Catch2)

add_executable(cppqedTests unit/testOperators.cc unit/testSlicing.cc unit/testHeaders.cc $<TARGET_OBJECTS:cppqedTestsHeadersCopy>)
target_link_libraries(cppqedTests PRIVATE CPPQED::CPPQEDelements Catch2::Catch2WithMain)

if(catch2_SOURCE_DIR)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
else()
  list(APPEND CMAKE_MODULE_PATH ${Catch2_DIR})
endif()
include(Catch)
catch_discover_tests(cppqedTests)

# Regression: jc_qjmc against the output of the code before the switch to cppqed-multidiagonal.
# Regenerate (only for an intended change of results): compareTrajectory.py <jc_qjmc> jc_qjmc.ref --regenerate
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME jc_qjmc_regression
           COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/regression/compareTrajectory.py
                   $<TARGET_FILE:jc_qjmc> ${CMAKE_CURRENT_SOURCE_DIR}/regression/jc_qjmc.ref)
else()
  message(WARNING "Python 3 not found: regression tests disabled")
endif()
