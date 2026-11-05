// Copyright András Vukics 2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \file
/// Stops with a message, before any other code of the framework runs, on a processor without the instruction set the
/// framework was built for (CMake option CPPQED_ARCH, #93). The code would otherwise die of an illegal instruction.
/// This file is compiled for the x86-64 baseline (CPPQEDutils/CMakeLists.txt), so that the check itself runs anywhere.
#include <cstdio>
#include <cstdlib>

#ifndef CPPQED_ARCH
#error "CPPQED_ARCH, the instruction-set level, is defined by CPPQEDutils/CMakeLists.txt; build this file through CMake"
#endif

namespace {

[[gnu::constructor(101)]] void checkInstructionSet()
{
  __builtin_cpu_init(); // this early, the processor may not have been detected yet
  if (!__builtin_cpu_supports(CPPQED_ARCH)) {
    std::fputs("C++QED was built for " CPPQED_ARCH " (CMake option CPPQED_ARCH), which this processor does not support.\n"
               "Rebuild for a lower level (-DCPPQED_ARCH=x86-64-v3 or x86-64-v2), or with -DCPPQED_ARCH= for the compiler's default.\n",
               stderr);
    std::_Exit(EXIT_FAILURE);
  }
}

} // anonymous
