cmake_minimum_required(VERSION 3.18)
project(ConditionsLib VERSION 0.1.0)

include_directories(include)

set(CMAKE_CXX_STANDARD 20)

# The macro trick is taken from:
# https://github.com/boostcon/cppnow_presentations_2017/blob/master/05-19-2017_friday/effective_cmake__daniel_pfeifer__cppnow_05-19-2017.pdf
set(as_subproject ConditionsLib)

macro(find_package)
  if (NOT "${ARGV0}" IN_LIST as_subproject)
    _find_package(${ARGV})
  endif()
endmacro()

add_subdirectory(src)
add_subdirectory(examples)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "Config.cmake.in"
  "ConditionsLibConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/
  PATH_VARS
    CMAKE_INSTALL_LIBDIR
  )

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/ConditionsLibConfigVersion.cmake
  VERSION 0.1.0
  COMPATIBILITY SameMajorVersion
  )

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/ConditionsLibConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/ConditionsLibConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ConditionsLib"
)

# Documentation
find_package(Doxygen)
if (DOXYGEN_FOUND)
  set(DOXYGEN_GENERATE_LATEX ON)
  doxygen_add_docs(
          doxygen
          ${PROJECT_SOURCE_DIR}
          COMMENT "Build Documentation"
  )
else()
  message("No doxygen documentation found.")
endif()