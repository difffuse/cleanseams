# SPDX-License-Identifier: GPL-3.0-or-later
include(CMakePackageConfigHelpers)
set_target_properties(cleanseams_core PROPERTIES EXPORT_NAME core)
set_target_properties(cleanseams_cli PROPERTIES EXPORT_NAME cli)
add_library(CleanSeams::core ALIAS cleanseams_core)
# A single runtime directory preserves the established extension/wheel loader
# contract on every platform. CMake imported targets locate the shared core here.
set(CS_INSTALL_PYTHONDIR "${CMAKE_INSTALL_BINDIR}")
foreach(directory CMAKE_INSTALL_BINDIR CMAKE_INSTALL_LIBDIR CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_DATADIR CS_INSTALL_PYTHONDIR)
  if(IS_ABSOLUTE "${${directory}}")
    message(FATAL_ERROR "${directory} must be relative for a relocatable CleanSeams SDK")
  endif()
endforeach()
if(APPLE)
  set(cs_origin "@loader_path")
  set_target_properties(cleanseams_core PROPERTIES INSTALL_NAME_DIR "@rpath")
elseif(UNIX)
  set(cs_origin "$ORIGIN")
endif()
if(DEFINED cs_origin)
  set_target_properties(cleanseams_core PROPERTIES INSTALL_RPATH "${cs_origin}")
  set_target_properties(cleanseams_cli PROPERTIES INSTALL_RPATH "${cs_origin}")
  if(CS_BUILD_BINDINGS)
    set_target_properties(cleanseams_native PROPERTIES INSTALL_RPATH "${cs_origin}")
  endif()
endif()
install(TARGETS cleanseams_core cleanseams_cli EXPORT CleanSeamsTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
if(CS_BUILD_BINDINGS)
  install(TARGETS cleanseams_native LIBRARY DESTINATION ${CS_INSTALL_PYTHONDIR}
    RUNTIME DESTINATION ${CS_INSTALL_PYTHONDIR})
endif()
install(DIRECTORY include/cleanseams DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.hpp")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/generated/cleanseams/export.hpp"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cleanseams)
# JSON is part of the public C++ API. Ship the same locked headers that built the
# core, without exporting a dependency on a build-tree target or system package.
install(DIRECTORY "${json_SOURCE_DIR}/include/nlohmann" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES COPYING deps.lock.json DESTINATION ${CMAKE_INSTALL_DATADIR}/cleanseams)
install(DIRECTORY LICENSES DESTINATION ${CMAKE_INSTALL_DATADIR}/cleanseams)
set(CS_RUNTIME_FILES "" CACHE STRING "Explicit app-local compiler runtimes")
set(CS_RUNTIME_LICENSE_FILES "" CACHE STRING "Notices for the selected compiler runtimes")
if(CS_RUNTIME_FILES AND NOT CS_RUNTIME_LICENSE_FILES)
  message(FATAL_ERROR "Provide CS_RUNTIME_LICENSE_FILES for the selected compiler runtimes")
endif()
foreach(runtime IN LISTS CS_RUNTIME_FILES CS_RUNTIME_LICENSE_FILES)
  if(NOT IS_ABSOLUTE "${runtime}" OR NOT EXISTS "${runtime}" OR IS_DIRECTORY "${runtime}")
    message(FATAL_ERROR "Runtime and notice inputs must be existing absolute files: ${runtime}")
  endif()
endforeach()
if(CS_RUNTIME_FILES)
  install(FILES ${CS_RUNTIME_FILES} DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
if(CS_RUNTIME_LICENSE_FILES)
  install(FILES ${CS_RUNTIME_LICENSE_FILES} DESTINATION ${CMAKE_INSTALL_DATADIR}/cleanseams/LICENSES/runtime)
endif()
configure_package_config_file(cmake/CleanSeamsConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/CleanSeamsConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CleanSeams
  PATH_VARS CS_INSTALL_PYTHONDIR)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/CleanSeamsConfigVersion.cmake"
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMinorVersion)
install(EXPORT CleanSeamsTargets NAMESPACE CleanSeams:: FILE CleanSeamsTargets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CleanSeams)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/CleanSeamsConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/CleanSeamsConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CleanSeams)
