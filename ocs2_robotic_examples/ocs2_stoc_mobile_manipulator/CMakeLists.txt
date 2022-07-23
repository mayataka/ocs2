cmake_minimum_required(VERSION 3.0.2)
project(ocs2_stoc_mobile_manipulator)

# Generate compile_commands.json for clang tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CATKIN_PACKAGE_DEPENDENCIES
  ocs2_core
  ocs2_oc
  ocs2_robotic_tools
  ocs2_mobile_manipulator
  ocs2_ipm
  ocs2_stoc
)

find_package(catkin REQUIRED COMPONENTS
  ${CATKIN_PACKAGE_DEPENDENCIES}
)

find_package(Boost REQUIRED COMPONENTS
  system
  filesystem
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

###################################
## catkin specific configuration ##
###################################

catkin_package(
  INCLUDE_DIRS
    include
    ${EIGEN3_INCLUDE_DIRS}
  LIBRARIES
    ${PROJECT_NAME}
  CATKIN_DEPENDS
    ${CATKIN_PACKAGE_DEPENDENCIES}
  DEPENDS
    Boost
)

###########
## Build ##
###########

# set(FLAGS
#   ${OCS2_CXX_FLAGS}
#   -Wno-ignored-attributes
#   -Wno-invalid-partial-specialization   # to silence warning with unsupported Eigen Tensor
# )

# Add directories for all targets
include_directories(
  include
  ${EIGEN3_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
)

# mobile maniulator interface library
add_library(${PROJECT_NAME}
  src/constraint/JointPositionLimitsIneq.cpp
  src/constraint/JointVelocityLimitsIneq.cpp
  src/STOC_MobileManipulatorInterface.cpp
)
add_dependencies(${PROJECT_NAME}
  ${catkin_EXPORTED_TARGETS}
)
target_link_libraries(${PROJECT_NAME}
  ${Boost_LIBRARIES}
  ${catkin_LIBRARIES}
)
target_compile_options(${PROJECT_NAME} PUBLIC ${FLAGS})

####################
## Clang tooling ###
####################

# find_package(cmake_clang_tools QUIET)
# if (cmake_clang_tools_FOUND)
#   message(STATUS "Run clang tooling")
#   add_clang_tooling(
#     TARGETS ${PROJECT_NAME}
#     SOURCE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/include
#     CT_HEADER_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
#     CF_WERROR
#   )
# endif (cmake_clang_tools_FOUND)

#############
## Install ##
#############

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)
install(DIRECTORY config
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
)

#############
## Testing ##
#############

# # Helper macro for adding target applications
# macro(add_ocs2_test APP_NAME APP_SRC)
#   catkin_add_gtest(${APP_NAME}
#     ${APP_SRC}
#   )
#   target_include_directories(${APP_NAME} PRIVATE
#     ${PROJECT_BINARY_DIR}/include
#   )
#   target_link_libraries(${APP_NAME}
#     gtest_main
#     ${PROJECT_NAME}
#     ${catkin_LIBRARIES}
#     ${Boost_LIBRARIES}
#   )
# endmacro()

