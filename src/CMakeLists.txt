add_library(ConditionsLib 
  FourierTable.cc
  SboxTools.cc
  VectorialBooleanFunction.cc
  basic_linear_algebra.cc
  enumeration_algorithm/BruteForce.cc
  enumeration_algorithm/BruteForceOptimizations.cc
  enumeration_algorithm/ExhaustiveBF.cc
  enumeration_algorithm/MinLeaves.cc)

include(GNUInstallDirs)

install(TARGETS ConditionsLib
  EXPORT ConditionsLibTargets
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include)
install(EXPORT ConditionsLibTargets
  FILE ConditionsLib.cmake
  NAMESPACE ConditionsLib::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ConditionsLib)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/include/ConditionsLib"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")


add_library(ConditionsLib::ConditionsLib ALIAS ConditionsLib)

