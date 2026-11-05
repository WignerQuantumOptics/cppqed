# Copyright András Vukics 2020–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring examples")

file(GLOB SCRIPTS . *.cc)
foreach(s ${SCRIPTS})
  get_filename_component(SCRIPT ${s} NAME_WE)
  add_executable(${SCRIPT} ${s})
  target_link_libraries(${SCRIPT} CPPQED::CPPQEDelements)
  if (TARGET PkgConfig::JEMALLOC)
    target_link_libraries(${SCRIPT} PkgConfig::JEMALLOC) # directly into the executable, so its malloc precedes the system's
  endif ()
endforeach(s)

