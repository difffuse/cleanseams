# Deterministic numerics also apply to compiled dependencies.
if(MSVC)
  add_compile_options(/fp:strict)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
else()
  add_compile_options(-fno-fast-math -ffp-contract=off)
endif()
# Keep checkout/build locations out of __FILE__ diagnostics and debug information.
# Directory options are inherited by dependencies, including sources fetched into
# the build tree. Preserve relative filenames so diagnostics remain useful.
# Native MSVC has no equivalent C/C++ -ffile-prefix-map option.
if(NOT MSVC)
  foreach(cs_language IN ITEMS C CXX)
    add_compile_options(
      "$<$<COMPILE_LANG_AND_ID:${cs_language},Clang,AppleClang,GNU>:-ffile-prefix-map=${CMAKE_SOURCE_DIR}=.>"
      "$<$<COMPILE_LANG_AND_ID:${cs_language},Clang,AppleClang,GNU>:-ffile-prefix-map=${CMAKE_BINARY_DIR}=.>")
  endforeach()
endif()
if(CS_SANITIZER STREQUAL "address")
  if(MSVC)
    add_compile_options(/fsanitize=address)
  else()
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined)
  endif()
elseif(CS_SANITIZER STREQUAL "thread" AND NOT MSVC)
  add_compile_options(-fsanitize=thread -fno-omit-frame-pointer)
  add_link_options(-fsanitize=thread)
elseif(NOT CS_SANITIZER STREQUAL "none")
  message(FATAL_ERROR "Unsupported sanitizer on this compiler: ${CS_SANITIZER}")
endif()
function(cs_project_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /permissive-)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
  endif()
endfunction()
