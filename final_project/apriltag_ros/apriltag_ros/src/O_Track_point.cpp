#include <iostream>
#include <ros/ros.h>
#include <Eigen/Core>
#include <mav_msgs/conversions.h>
#include <apriltag_ros/AprilTagDetectionArray.h>
#include <trajectory_msgs/MultiDOFJointTrajectory.h>

using namespace std;

trajectory_msgs::MultiDOFJointTrajectory HoldPoint;

void HoldCallback(const apriltag_ros::AprilTagDetectionArray::ConstPtr &msg)
{
  HoldPoint = msg;
}


int main(int argc, char **argv)
{

  // Node initialization
  ros::init(argc, argv, "Tracker");
  ros::NodeHandle n;

  multirotros_to_jusky NC = multirotros_to_jusky(&n);
  NC.Get_apriltag_pose();

  boost::shared_ptr<apriltag_ros::AprilTagDetectionArray const> sharedEdge;
  ros::Publisher trajectory_husky_pub;

  // Node Assignment
  trajectory_husky_pub = n.advertise<trajectory_msgs::MultiDOFJointTrajectory>("/Husky_point", 1);

  sharedEdge = ros::topic::waitForMessage<apriltag_ros::AprilTagDetectionArray>("/tag_detections", n);

  apriltag_ros::AprilTagDetectionArray Husky_pose;

  // position msg
  if(sharedEdge != NULL)
  {
    ROS_INFO("True");
    Husky_pose = *sharedEdge;
  }
  
  trajectory_msgs::MultiDOFJointTrajectory trajectory_husky_msg;

  // Desired position
  float x,y,z;
  x = Husky_pose.detections[0].pose.pose.pose.position.x;
  y = Husky_pose.detections[0].pose.pose.pose.position.y;
  z = Husky_pose.detections[0].pose.pose.pose.position.z;

  float Real_x, Real_y, Real_z;
  Real_x = y + 15.5;
  Real_y = -(x+2) + 10;
  Real_z = z;

  Eigen::Vector3d desired_position(Real_x, Real_y, Real_z);

  // Desired attitude
  const float DEG_2_RAD = M_PI / 180.0;

  double desired_yaw = 2 * DEG_2_RAD;

  mav_msgs::msgMultiDofJointTrajectoryFromPositionYaw(desired_position, desired_yaw, &trajectory_husky_msg);

  trajectory_husky_pub.publish(trajectory_husky_msg);

  return 0;
}