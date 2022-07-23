cmake_minimum_required(VERSION 3.0.2)
project(ocs2_stoc_mobile_manipulator_ros)

# Generate compile_commands.json for clang tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CATKIN_PACKAGE_DEPENDENCIES
  roslib
  tf
  urdf
  kdl_parser
  robot_state_publisher
  visualization_msgs
  geometry_msgs
  ocs2_ros_interfaces
  ocs2_core
  ocs2_stoc
  ocs2_mpc
  ocs2_robotic_tools
  ocs2_robotic_assets
  ocs2_pinocchio_interface
  ocs2_self_collision
  ocs2_self_collision_visualization
  ocs2_mobile_manipulator
  ocs2_mobile_manipulator_ros
  ocs2_stoc_mobile_manipulator
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
    ${EIGEN3_INCLUDE_DIRS}
  CATKIN_DEPENDS
    ${CATKIN_PACKAGE_DEPENDENCIES}
  DEPENDS
    Boost
)

###########
## Build ##
###########

set(FLAGS
  ${OCS2_CXX_FLAGS}
  -Wno-ignored-attributes
  -Wno-invalid-partial-specialization   # to silence warning with unsupported Eigen Tensor
)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
)

# MPC node
add_executable(mobile_manipulator_mpc_node
  src/MobileManipulatorMpcNode.cpp
)
add_dependencies(mobile_manipulator_mpc_node
  ${catkin_EXPORTED_TARGETS}
)
target_link_libraries(mobile_manipulator_mpc_node
  ${catkin_LIBRARIES}
)
target_compile_options(mobile_manipulator_mpc_node PUBLIC ${FLAGS})

####################
## Clang tooling ###
####################

# find_package(cmake_clang_tools QUIET)
# if (cmake_clang_tools_FOUND)
#   message(STATUS "Run clang tooling")
#   add_clang_tooling(
#     TARGETS mobile_manipulator_mpc_node mobile_manipulator_dummy_mrt_node
#     SOURCE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/include
#     CT_HEADER_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
#     CF_WERROR
#   )
# endif (cmake_clang_tools_FOUND)

#############
## Install ##
#############

install(
  TARGETS 
    mobile_manipulator_mpc_node 
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
install(DIRECTORY launch 
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
)
