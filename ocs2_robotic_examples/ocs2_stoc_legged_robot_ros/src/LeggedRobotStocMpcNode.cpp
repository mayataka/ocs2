#include <ros/init.h>
#include <ros/package.h>

#include <ocs2_legged_robot/LeggedRobotInterface.h>
#include <ocs2_stoc_legged_robot/STOC_LeggedRobotInterface.h>
#include <ocs2_ros_interfaces/mpc/MPC_ROS_Interface.h>
#include <ocs2_ros_interfaces/synchronized_module/RosReferenceManager.h>
#include <ocs2_stoc/STOC_MPC.h>

#include "ocs2_legged_robot_ros/gait/GaitReceiver.h"

using namespace ocs2;
using namespace legged_robot;

int main(int argc, char** argv) {
  const std::string robotName = "legged_robot";

  // Initialize ros node
  ros::init(argc, argv, robotName + "_mpc");
  ros::NodeHandle nodeHandle;
  // Get node parameters
  std::string stocFile, taskFile, urdfFile, referenceFile;
  nodeHandle.getParam("/stocFile", stocFile);
  nodeHandle.getParam("/taskFile", taskFile);
  nodeHandle.getParam("/urdfFile", urdfFile);
  nodeHandle.getParam("/referenceFile", referenceFile);

  // Robot interface
  // LeggedRobotInterface interface(taskFile, urdfFile, referenceFile);
  STOC_LeggedRobotInterface stocInterface(stocFile, taskFile, urdfFile, referenceFile);

  // Gait receiver
  auto gaitReceiverPtr =
      std::make_shared<GaitReceiver>(nodeHandle, stocInterface.getSwitchedModelReferenceManagerPtr()->getGaitSchedule(), robotName);

  // ROS ReferenceManager
  auto rosReferenceManagerPtr = std::make_shared<RosReferenceManager>(robotName, stocInterface.getReferenceManagerPtr());
  rosReferenceManagerPtr->subscribe(nodeHandle);

  // MPC
  STOC_MPC mpc(stocInterface.mpcSettings(), stocInterface.stocSettings(), stocInterface.getOptimalControlProblem(),
               stocInterface.getInitializer());
  mpc.getSolverPtr()->setReferenceManager(rosReferenceManagerPtr);
  mpc.getSolverPtr()->addSynchronizedModule(gaitReceiverPtr);

  // Launch MPC ROS node
  MPC_ROS_Interface mpcNode(mpc, robotName);
  mpcNode.launchNodes(nodeHandle);

  // Successful exit
  return 0;
}
