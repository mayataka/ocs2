cmake_minimum_required(VERSION 3.0.2)
project(ocs2_stoc)

set(CATKIN_PACKAGE_DEPENDENCIES
  ocs2_core
  ocs2_oc
  ocs2_sqp
  ocs2_sto
  ocs2_ipm
  ocs2_ipm_oc
  ocs2_mpc
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
  CATKIN_DEPENDS
    ${CATKIN_PACKAGE_DEPENDENCIES}
  LIBRARIES
    ${PROJECT_NAME}
  DEPENDS
    Boost
)

###########
## Build ##
###########

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
)

# ocs2_stoc library
add_library(${PROJECT_NAME}
  src/sto/LinearQuadraticApproximator.cpp
  src/sto/IpmVariables.cpp
  src/sto/IpmVariablesEliminator.cpp
  src/sto/IpmVariablesRetriver.cpp
  src/riccati_recursion/BackwardRiccatiRecursion.cpp
  src/riccati_recursion/ForwardRiccatiRecursion.cpp
  src/riccati_recursion/RiccatiRecursion.cpp
  src/TimeDiscretization.cpp
  src/STOC_Settings.cpp
  src/STOC.cpp
)
add_dependencies(${PROJECT_NAME}
  ${catkin_EXPORTED_TARGETS}
)
target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
)
target_compile_options(${PROJECT_NAME} PUBLIC ${OCS2_CXX_FLAGS})

#########################
###   CLANG TOOLING   ###
#########################
find_package(cmake_clang_tools QUIET)
if(cmake_clang_tools_FOUND)
  message(STATUS "Run clang tooling for target " ${PROJECT_NAME})
  add_clang_tooling(
    TARGETS
    ${PROJECT_NAME}
    SOURCE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/include
    CT_HEADER_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
    CF_WERROR
  )
endif(cmake_clang_tools_FOUND)

#############
## Install ##
#############
install(
  TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

install(DIRECTORY include/${PROJECT_NAME}/
        DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION})

#############
## Testing ##
#############

catkin_add_gtest(test_${PROJECT_NAME}
  test/testTimeDiscretization.cpp
)
add_dependencies(test_${PROJECT_NAME} ${catkin_EXPORTED_TARGETS})
target_link_libraries(test_${PROJECT_NAME}
  ${PROJECT_NAME}
  ${catkin_LIBRARIES}
  gtest_main
)

catkin_add_gtest(Exp0Test_${PROJECT_NAME}
  test/Exp0Test.cpp
)
add_dependencies(Exp0Test_${PROJECT_NAME} ${catkin_EXPORTED_TARGETS})
target_link_libraries(Exp0Test_${PROJECT_NAME} 
  ${PROJECT_NAME}
  ${catkin_LIBRARIES}
  gtest_main
)

catkin_add_gtest(Exp1Test_${PROJECT_NAME}
  test/Exp1Test.cpp
)
add_dependencies(Exp1Test_${PROJECT_NAME} ${catkin_EXPORTED_TARGETS})
target_link_libraries(Exp1Test_${PROJECT_NAME} 
  ${PROJECT_NAME}
  ${catkin_LIBRARIES}
  gtest_main
)

catkin_add_gtest(testCircularKinematics_${PROJECT_NAME}
  test/testCircularKinematics.cpp
)
add_dependencies(testCircularKinematics_${PROJECT_NAME} ${catkin_EXPORTED_TARGETS})
target_link_libraries(testCircularKinematics_${PROJECT_NAME} 
  ${PROJECT_NAME}
  ${catkin_LIBRARIES}
  gtest_main
)