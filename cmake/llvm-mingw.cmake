# SPDX-License-Identifier: GPL-3.0-or-later
# Cross-compile Windows x86-64 with an extracted LLVM-MinGW toolchain.
# Pass -DLLVM_MINGW=/absolute/toolchain/path.
set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_SYSTEM_PROCESSOR AMD64)
set(LLVM_MINGW "" CACHE PATH "Extracted LLVM-MinGW distribution")
if(NOT IS_DIRECTORY "${LLVM_MINGW}/x86_64-w64-mingw32")
  message(FATAL_ERROR "Pass -DLLVM_MINGW pointing to the extracted toolchain")
endif()
list(APPEND CMAKE_TRY_COMPILE_PLATFORM_VARIABLES LLVM_MINGW)
set(CMAKE_C_COMPILER "${LLVM_MINGW}/bin/x86_64-w64-mingw32-clang")
set(CMAKE_CXX_COMPILER "${LLVM_MINGW}/bin/x86_64-w64-mingw32-clang++")
set(CMAKE_RC_COMPILER "${LLVM_MINGW}/bin/x86_64-w64-mingw32-windres")
set(CMAKE_FIND_ROOT_PATH "${LLVM_MINGW}/x86_64-w64-mingw32")
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)
