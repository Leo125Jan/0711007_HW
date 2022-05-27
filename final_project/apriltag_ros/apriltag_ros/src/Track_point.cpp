#include <iostream>
#include <unistd.h>
#include <ros/ros.h>
#include <Eigen/Core>
#include <std_msgs/Int64.h>
#include <mav_msgs/conversions.h>
#include <apriltag_ros/AprilTagDetectionArray.h>
#include <trajectory_msgs/MultiDOFJointTrajectory.h>

using namespace std;

class multirotros_to_jusky
{
  private:

    // Node Initialization
    ros::Publisher trajectory_husky_pub_d;
    ros::Publisher trajectory_husky_pub_1;
    boost::shared_ptr<apriltag_ros::AprilTagDetectionArray const> sharedEdge;

    apriltag_ros::AprilTagDetectionArray Husky_pose;
    trajectory_msgs::MultiDOFJointTrajectory trajectory_husky_msg_d;
    trajectory_msgs::MultiDOFJointTrajectory trajectory_husky_msg_1;

    float x, y, z, Real_x_d, Real_y_d, Real_z_d, Real_x_1, Real_y_1, Real_z_1;
    Eigen::Vector3d desired_position_d;
    Eigen::Vector3d desired_position_1;

    const float DEG_2_RAD = M_PI / 180.0;
    double desired_yaw = 2 * DEG_2_RAD;

  public:

    multirotros_to_jusky(ros::NodeHandle *n)
    {
      trajectory_husky_pub_d = n->advertise<trajectory_msgs::MultiDOFJointTrajectory>("/Husky_point", 1);
      trajectory_husky_pub_1 = n->advertise<trajectory_msgs::MultiDOFJointTrajectory>("/Husky_point_to_1", 1);  
    }

    void Get_apriltag_pose_d(void)
    {
      sharedEdge = ros::topic::waitForMessage<apriltag_ros::AprilTagDetectionArray>("/tag_detections");

      if(sharedEdge != NULL)
      {
        Husky_pose = *sharedEdge;
      }

      // default
      x = Husky_pose.detections[0].pose.pose.pose.position.x;
      y = Husky_pose.detections[0].pose.pose.pose.position.y;
      z = Husky_pose.detections[0].pose.pose.pose.position.z;

      Real_x_d = -(y + 0.36) + 15.5;
      Real_y_d = -(x+2) + 10;
      Real_z_d = z - 10;

      desired_position_d << Real_x_d, Real_y_d, Real_z_d;

      mav_msgs::msgMultiDofJointTrajectoryFromPositionYaw(desired_position_d, desired_yaw, &trajectory_husky_msg_d);

      trajectory_husky_pub_d.publish(trajectory_husky_msg_d);

      // One
      Real_x_1 = y + 15.5;
      Real_y_1 = -(x+2) + 10;
      Real_z_1 = z - 5;

      desired_position_1 << Real_x_1, Real_y_1, Real_z_1;

      mav_msgs::msgMultiDofJointTrajectoryFromPositionYaw(desired_position_1, desired_yaw, &trajectory_husky_msg_1);

      trajectory_husky_pub_1.publish(trajectory_husky_msg_1);
    }
};

int main(int argc, char **argv)
{
  ros::init(argc, argv, "Tracker");
  ros::NodeHandle n;

  multirotros_to_jusky NC = multirotros_to_jusky(&n);
  NC.Get_apriltag_pose_d();

  return 0;
}

  // <!-- Robot two-->
  // <group ns="firefly_1">

  //   <include file="$(find rotors_gazebo)/launch/spawn_mav.launch">

  //     <arg name="mav_name" value="$(arg mav_name)" />
  //     <arg name="namespace" value="firefly_1"/>
  //     <arg name="model" value="$(find rotors_description)/urdf/mav_with_vi_sensor.gazebo" />
  //     <arg name="enable_logging" value="$(arg enable_logging)" />
  //     <arg name="enable_ground_truth" value="$(arg enable_ground_truth)" />
  //     <arg name="log_file" value="$(arg log_file)"/>
  //     <arg name="y" value="10.0"/>
  //     <arg name="x" value="10.0"/>

  //   </include>
  // </group>

    //   <group ns="firefly_1">

    //     <node name="lee_position_controller_node" pkg="rotors_control" type="lee_position_controller_node" output="screen">
    //         <rosparam command="load" file="$(find rotors_gazebo)/resource/lee_controller_$(arg mav_name).yaml" />
    //         <rosparam command="load" file="$(find rotors_gazebo)/resource/$(arg mav_name).yaml" />
    //         <remap from="odometry" to="odometry_sensor1/odometry" />
    //     </node>
    //     <node name="waypoint_publisher_1" pkg="rotors_gazebo" type="waypoint_publisher_1" output="screen" args="15 -10 10 0 2"/>
    //     <node name="robot_state_publisher" pkg="robot_state_publisher" type="robot_state_publisher" />
    //     <node name="joint_state_publisher" pkg="joint_state_publisher" type="joint_state_publisher" />
        
    // </group>