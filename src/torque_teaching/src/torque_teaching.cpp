
#include "ros/ros.h"
#include "std_msgs/Float64.h"
#include <kdl_parser/kdl_parser.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <kdl/chain.hpp>
#include <kdl/chainfksolver.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/frames_io.hpp>
#include <stdio.h>
#include <kdl/kdl.hpp>
#include <kdl/jntarray.hpp>
#include <sensor_msgs/JointState.h>
#include <angles/angles.h>
#include <Eigen/Dense>
#include <interactive_markers/interactive_marker_server.h>
#include <boost/make_shared.hpp>
#include <tf/tf.h>
#include <tf_conversions/tf_kdl.h>
#include <kdl/chainiksolverpos_nr.hpp>

#include <math.h>
#include "ikbasedonqp.h"
#include <SPLINTER/bspline.h>
#include <SPLINTER/datatable.h>
#include <SPLINTER/bsplinebasis.h>
#include "trajectory_generate.h"
#include <kdl/joint.hpp>
#include <robotiq_force_torque_sensor/ft_sensor.h>
#include <robotiq_force_torque_sensor/sensor_accessor.h>

using namespace KDL;


#define PI 3.14159

KDL::JntArray jointpositions_now(6);
KDL::JntArray jointpositions_now_first(6);
KDL::JntArray jointpositions_gazebo(6);
KDL::JntArray jointpositions(6);

unsigned int j,m=0;
double alfa_temp;
double beta_temp;
double gammar_temp;
KDL::Frame cartpos_temp;
bool isGetRealJointValue = false;
typedef Eigen::Matrix<double,6,6> Matrix6d;
typedef Eigen::Matrix<double,6,1> Matrix61d;
typedef Eigen::Matrix<double,3,3> Matrix33d;
typedef Eigen::Matrix<double,3,1> Matrix31d;
double servoj_time = 0.003;
double t_homing = 5;
double t=0;
std::vector<float> force_torque(6,0.0);

std::vector<std::vector<float> > force_torque_calibration(9,force_torque);
Matrix31d L;
KDL::Chain chain_robotiq;
Matrix61d des_speed;
Matrix61d des_position;

Matrix61d delta_position;

unsigned int computeDepthFromRoot(KDL::SegmentMap::const_iterator el, const KDL::SegmentMap::const_iterator &root)
{
    unsigned int result = 0;
    while (el != root) {
    ++result;
    el = el->second.parent;
    }
    return result;
}

const std::string& findEndeffector(const KDL::Tree &tree)
{
    const KDL::SegmentMap &segments = tree.getSegments();
    unsigned int maxDepth = 0;
    KDL::SegmentMap::const_iterator eef = segments.end();

        for (KDL::SegmentMap::const_iterator it = segments.begin(), end = segments.end(); it != end; ++it)
        {
            if (it->second.children.size() == 0)
            {
                unsigned int depth = computeDepthFromRoot(it, tree.getRootSegment());

           if (depth > maxDepth || eef == segments.end())
           {
                eef = it;
                maxDepth = depth;
           }
            }
        }
    return eef->first;
}

/*
  return the jointpositions at t time after the cubic interpolation
 */
KDL::JntArray homing(KDL::JntArray jointpositions_s,double t,double t_s,double t_f,
                     std::vector<ros::Publisher> chatter_pub)
{
    if(t<=t_homing){

        KDL::JntArray jointpositions_standard(6);
        KDL::JntArray jointpositions_f(6);
//        jointpositions_standard.data <<0,-PI/2,PI/2,0,0,0;
        jointpositions_standard.data <<-(91.71/180)*PI,(-98.96/180)*PI,(-126.22/180)*PI,(-46.29/180)*PI,(91.39/180)*PI,(-1.78/180)*PI;

        std_msgs::Float64 msg1;
        std_msgs::Float64 msg2;
        std_msgs::Float64 msg3;
        std_msgs::Float64 msg4;
        std_msgs::Float64 msg5;
        std_msgs::Float64 msg6;
        sensor_msgs::JointState jointpositions_msg;
        jointpositions_f = jointpositions_standard;

        for(unsigned int i=0;i<6;i++){

            jointpositions(i) = jointpositions_s(i)
             +(3/((t_f-t_s)*(t_f-t_s)))*(jointpositions_f(i)-jointpositions_s(i))*(t-t_s)*(t-t_s)
             -(2/((t_f-t_s)*(t_f-t_s)*(t_f-t_s)))*(jointpositions_f(i)-jointpositions_s(i))*(t-t_s)*(t-t_s)*(t-t_s);
        }


        msg1.data=jointpositions(0);
        msg2.data=jointpositions(1);
        msg3.data=jointpositions(2);
        msg4.data=jointpositions(3);
        msg5.data=jointpositions(4);
        msg6.data=jointpositions(5);

            for(unsigned int i=0;i<6;i++){
                jointpositions_msg.position.push_back(jointpositions(i));
            }
        //to gazebo
        chatter_pub[0].publish(msg1);
        chatter_pub[1].publish(msg2);
        chatter_pub[2].publish(msg3);
        chatter_pub[3].publish(msg4);
        chatter_pub[4].publish(msg5);
        chatter_pub[5].publish(msg6);

        //to real robot
        chatter_pub[6].publish(jointpositions_msg);
    }
    return jointpositions;
}


void realrobotcallback_JointStates(const sensor_msgs::JointState::ConstPtr& realrobot_jointstates){

    j++;
        for(unsigned int i=0;i<6;i++){
            jointpositions_now(i) = realrobot_jointstates->position[i];
        }
        if(j==1){
          jointpositions_now_first = jointpositions_now;
        }
//    std::cout<<"ok"<<std::endl;
    isGetRealJointValue = true;

}

void gazeborobotcallback_JointStates(const sensor_msgs::JointState::ConstPtr& gazeborobot_jointstates){

    for(unsigned int i=0;i<6;i++){
        jointpositions_gazebo(i) = gazeborobot_jointstates->position[i];
    }

}



//t_user is the user time.t is the system time

KDL::JntArray trajectory_do(double t,std::vector<double> T_user,std::vector<std::vector<double> > cartpospoint,
                   KDL::Chain my_chain,std::vector<ros::Publisher> chatter_pub){

    std::vector<double> T;
    int n= T_user.size();

    for(unsigned int i=0;i<n;i++){
        T.push_back(T_user[i]+t_homing);
    }

    if(t>=T.front()&&t<=T.back()){

        std_msgs::Float64 msg1;
        std_msgs::Float64 msg2;
        std_msgs::Float64 msg3;
        std_msgs::Float64 msg4;
        std_msgs::Float64 msg5;
        std_msgs::Float64 msg6;
        sensor_msgs::JointState jointpositions_msg;

        if(T.size()>2){

            std::cout<<"ok1111"<<std::endl;

            trajectory_generate tg;
            tg.spline_generate(T,cartpospoint);

            Matrix61d bound_low;
            Matrix61d bound_upper;
            Matrix61d des_cartpos;
            Matrix61d des_cartendvel;

            des_cartpos << tg.s_x(t),tg.s_y(t),tg.s_z(t),
                        tg.s_gammar(t),tg.s_beta(t),tg.s_alfa(t);

            des_cartendvel << tg.s_x.deriv(1,t),tg.s_y.deriv(1,t),tg.s_z.deriv(1,t),
                              tg.s_gammar.deriv(1,t),tg.s_beta.deriv(1,t),tg.s_alfa.deriv(1,t);


            bound_low << -2*PI,-2*PI,-2*PI,-2*PI,-2*PI,-2*PI;
            bound_upper << 2*PI,2*PI ,2*PI,2*PI,2*PI,2*PI;

        // IK

            jointpositions = IKbasedonQP::IK(jointpositions,jointpositions_now,des_cartpos,des_cartendvel,
                                   my_chain,bound_low,bound_upper);


          }


        else{


            std::vector<double> T_new;
            std::vector<std::vector<double>> cartpospoint_new;
            double T_mid = (T[0]+T[1])/2;

            std::vector<double> cartpos_mid(6,0);
            cartpos_mid[0]=(cartpospoint[0][0]+cartpospoint[1][0])/2;
            cartpos_mid[1]=(cartpospoint[0][1]+cartpospoint[1][1])/2;
            cartpos_mid[2]=(cartpospoint[0][2]+cartpospoint[1][2])/2;
            cartpos_mid[3]=(cartpospoint[0][3]+cartpospoint[1][3])/2;
            cartpos_mid[4]=(cartpospoint[0][4]+cartpospoint[1][4])/2;
            cartpos_mid[5]=(cartpospoint[0][5]+cartpospoint[1][5])/2;

            T_new.push_back(T[0]);
            T_new.push_back(T_mid);
            T_new.push_back(T[1]);

            cartpospoint_new.push_back(cartpospoint[0]);
            cartpospoint_new.push_back(cartpos_mid);
            cartpospoint_new.push_back(cartpospoint[1]);

            trajectory_generate tg;
            tg.spline_generate(T_new,cartpospoint_new);

            Matrix61d bound_low;
            Matrix61d bound_upper;
            Matrix61d des_cartpos;
            Matrix61d des_cartendvel;

            des_cartpos << tg.s_x(t),tg.s_y(t),tg.s_z(t),
                        tg.s_gammar(t),tg.s_beta(t),tg.s_alfa(t);

            des_cartendvel << tg.s_x.deriv(1,t),tg.s_y.deriv(1,t),tg.s_z.deriv(1,t),
                              tg.s_gammar.deriv(1,t),tg.s_beta.deriv(1,t),tg.s_alfa.deriv(1,t);


            bound_low << -2*PI,-2*PI,-2*PI,-2*PI,-2*PI,-2*PI;
            bound_upper << 2*PI,2*PI ,2*PI,2*PI,2*PI,2*PI;

        // IK

            jointpositions = IKbasedonQP::IK(jointpositions,jointpositions_now,des_cartpos,des_cartendvel,
                                   my_chain,bound_low,bound_upper);


        }


        msg1.data=jointpositions(0);
        msg2.data=jointpositions(1);
        msg3.data=jointpositions(2);
        msg4.data=jointpositions(3);
        msg5.data=jointpositions(4);
        msg6.data=jointpositions(5);

        for(unsigned int i=0;i<6;i++){
          jointpositions_msg.position.push_back(jointpositions(i));
        }

        //to gazebo
        chatter_pub[0].publish(msg1);
        chatter_pub[1].publish(msg2);
        chatter_pub[2].publish(msg3);
        chatter_pub[3].publish(msg4);
        chatter_pub[4].publish(msg5);
        chatter_pub[5].publish(msg6);

        //to real robot
        chatter_pub[6].publish(jointpositions_msg);

        // ROSINFO
        std::cout<<"jointposition_real "
                << jointpositions_now(0)*180/PI<<" "
                << jointpositions_now(1)*180/PI<<" "
                << jointpositions_now(2)*180/PI<<" "
                << jointpositions_now(3)*180/PI<<" "
                << jointpositions_now(4)*180/PI<<" "
                << jointpositions_now(5)*180/PI<<" "
                << std::endl;
        std::cout<<"jointposition_plan "
                << jointpositions(0)*180/PI<<" "
                << jointpositions(1)*180/PI<<" "
                << jointpositions(2)*180/PI<<" "
                << jointpositions(3)*180/PI<<" "
                << jointpositions(4)*180/PI<<" "
                << jointpositions(5)*180/PI<<" "
                << std::endl;
    }


    return jointpositions;

}


void trajectory_sleep(double t,double t_now,double sleeptime){

    double t_start;
    double t_end;
    t_start= t_now;
    t_end = t_start+sleeptime;
    if(t>=t_start&&t<=t_end){

         ROS_INFO("sleep");

       }
}

void reCallback(const robotiq_force_torque_sensor::ft_sensor& msg)
{
    ROS_INFO("I heard: FX[%f] FY[%f] FZ[%f] MX[%f] MY[%f] MZ[%f]", msg.Fx,msg.Fy,msg.Fz,msg.Mx,msg.My,msg.Mz);
    force_torque[0]=msg.Fx;
    force_torque[1]=msg.Fy;
    force_torque[2]=msg.Fz;
    force_torque[3]=msg.Mx;
    force_torque[4]=msg.My;
    force_torque[5]=msg.Mz;
    ROS_INFO("I heard: force_torque", force_torque[0],force_torque[1],force_torque[2],force_torque[3],
            force_torque[4],force_torque[5]);

}


/**
 * \brief this function is used to calibration the center of the end effector
 * @param force_torque_calibration - the data from the torque sensor
   */
void centroid_weight_calibration(std::vector<std::vector<float> > force_torque_calibration)
{

    Matrix33d F_init;
    Matrix31d M_init;
    F_init<<0,0,0,0,0,0,0,0,0;
    M_init<<0,0,0;
    Matrix33d sum_FF;
    Matrix31d sum_FM;
    Matrix33d FF_inverse;
    std::vector<Matrix33d> F_trans(9,F_init);
    std::vector<Matrix33d> F(9,F_init);
    std::vector<Matrix31d> M(9,M_init);

    //debug
    std::cout<<"force_torque_calibration"<<std::endl;
    std::cout<<force_torque_calibration[0][0]<<std::endl;
    std::cout<<force_torque_calibration[1][0]<<std::endl;
    std::cout<<force_torque_calibration[2][0]<<std::endl;
    std::cout<<force_torque_calibration[3][0]<<std::endl;
    std::cout<<force_torque_calibration[4][0]<<std::endl;
    std::cout<<force_torque_calibration[5][0]<<std::endl;
    std::cout<<force_torque_calibration[6][0]<<std::endl;
    std::cout<<force_torque_calibration[7][0]<<std::endl;
    std::cout<<force_torque_calibration[8][0]<<std::endl;

    for(unsigned i=0;i<9;i++){
        F[i]<<0,-force_torque_calibration[i][2],force_torque_calibration[i][1],
               force_torque_calibration[i][2],0,-force_torque_calibration[i][0],
              -force_torque_calibration[i][1],force_torque_calibration[i][0],0;

    }

    //debug
    std::cout<<"F[i]:"<<std::endl;
    std::cout<<F[0]<<std::endl;
    std::cout<<F[1]<<std::endl;
    std::cout<<F[2]<<std::endl;

    for(unsigned i=0;i<9;i++){
        M[i]<<force_torque_calibration[i][3],force_torque_calibration[i][4],force_torque_calibration[i][5];
    }

    for(unsigned i=0;i<9;i++){
        F_trans[i]=F[i].transpose();
    }
    //debug
    std::cout<<"F_trans:"<<std::endl;
    std::cout<<F_trans[0]<<std::endl;
    std::cout<<F_trans[1]<<std::endl;
    std::cout<<F_trans[2]<<std::endl;

    sum_FF = (F_trans[0]*F[0])+(F_trans[1]*F[1])+(F_trans[2]*F[2])
              +(F_trans[3]*F[3])+(F_trans[4]*F[4])+(F_trans[5]*F[5])
              +(F_trans[6]*F[6])+(F_trans[7]*F[7])+(F_trans[8]*F[8]);
    //debug
    std::cout<<"sum_FF:"<<std::endl;
    std::cout<<sum_FF<<std::endl;

    FF_inverse = sum_FF.inverse();
    //debug
    std::cout<<"FF_inverse:"<<std::endl;
    std::cout<<FF_inverse<<std::endl;

    sum_FM = F_trans[0]*M[0]+F_trans[1]*M[1]+F_trans[2]*M[2]
            +F_trans[3]*M[3]+F_trans[4]*M[4]+F_trans[5]*M[5]
            +F_trans[6]*M[6]+F_trans[7]*M[7]+F_trans[8]*M[8];
    //debug
    std::cout<<"sum_FM:"<<std::endl;
    std::cout<<sum_FM<<std::endl;

    L = FF_inverse*sum_FM;

    //debug
    std::cout<<"L "<<L<<std::endl;
}

/**
 * \brief this function is used to do the transfom between force and postion speed
 * @param F - force after compensatation
   */
void mapping(Matrix61d F){
    Matrix61d low_limit;
    Matrix61d up_limit;

    double k1;
    double k2;

    //debug k1 k2 low_limit up_limit
    k1 = 0.1;
    k2 = 0.1;
    low_limit<<0.5,0.5,0.5,0.5,0.5,0.5;
    up_limit<<10,10,10,10,10,10;

    for(unsigned i=0;i<6;i++){
        if(fabs(F[i])>10&&fabs(F[i])<20){
            des_speed[i] = k1*F[i];
            delta_position[i] =k2*F[i];
        }
        else{
            des_speed[i] = 0;
            delta_position[i] =0;
        }
    }

}

/**
 * \brief it can realize the force-chasing function
 * @param jointpositions - position that will send to the robot
 * @param my_chain - chain of ur
 * @param chatter_pub - ros publisher
   */
KDL::JntArray force_chasing(KDL::JntArray jointpositions,KDL::Chain my_chain,
                            std::vector<ros::Publisher> chatter_pub){

    //ros message
    sensor_msgs::JointState jointpositions_msg;
    Matrix61d bound_low;
    Matrix61d bound_upper;
    bound_low << -2*PI,-2*PI,-2*PI,-2*PI,-2*PI,-2*PI;
    bound_upper << 2*PI,2*PI ,2*PI,2*PI,2*PI,2*PI;

//    // caculate the Jacobian of robotiq torque sensor
//    KDL::ChainJntToJacSolver ur_JacSlover = KDL::ChainJntToJacSolver(chain_robotiq);
//    KDL::JntArray joint_in(6);
//    KDL::Jacobian ur_robotiq_jac(6);
//    joint_in = jointpositions_now;
//    ur_JacSlover.JntToJac(joint_in, ur_robotiq_jac, -1);
//    std::cout<<"ur_robotiq_jac_Jacobian "<< std::endl << ur_robotiq_jac.data << std::endl;

    // force compensation
    Matrix61d F0;
    Matrix61d F0_after;
    Matrix61d F_end;
    Matrix6d H;
    double G;
//    Matrix6d j_robotiq;


//    j_robotiq = ur_robotiq_jac.data;
    F_end<< force_torque[0],force_torque[1],force_torque[2],
            force_torque[3],force_torque[4],force_torque[5];
    std::cout<<"F_end "<<F_end<<std::endl;
    G = sqrt(force_torque[0]*force_torque[0]+force_torque[1]*force_torque[1]+force_torque[2]*force_torque[2]);



    //FK
    Matrix61d now_position;
    KDL::Frame frame_ur;
    KDL::Frame ur_robotiq;
    KDL::Frame frame_robotiq;
    KDL::Frame frame_L;
    KDL::Frame robotiq_L;


    double alfa_robotiq;
    double beta_robotiq;
    double gammar_robotiq;
    double alfa_ur;
    double beta_ur;
    double gammar_ur;

    frame_ur = IKbasedonQP::FK(jointpositions_now,my_chain);
    frame_ur.M.GetEulerZYX(alfa_ur,beta_ur,gammar_ur);
    ur_robotiq = Frame(Rotation::RPY(0,PI/2,PI),Vector(0.03,0,0));
    robotiq_L = Frame(Rotation::RPY(0,0,0),Vector(L[0],L[1],L[2]));
    frame_robotiq = frame_ur*ur_robotiq;
    frame_L = frame_robotiq*robotiq_L;
    frame_robotiq.M.GetEulerZYX(alfa_robotiq,beta_robotiq,gammar_robotiq);
    Matrix33d R;

    R<< frame_robotiq.M.data[0],frame_robotiq.M.data[1],frame_robotiq.M.data[2],
        frame_robotiq.M.data[3],frame_robotiq.M.data[4],frame_robotiq.M.data[5],
        frame_robotiq.M.data[6],frame_robotiq.M.data[7],frame_robotiq.M.data[8];

    Matrix33d S_p;

    S_p<<0,-frame_robotiq.p(2),frame_robotiq.p(1),
         frame_robotiq.p(2),0,-frame_robotiq.p(0),
         -frame_robotiq.p(1),frame_robotiq.p(0),0;

    H.block<3,3>(0,0) = R.transpose();
    H.block<3,3>(0,3) = -(R.transpose())*S_p;
    H.block<3,3>(3,0) = Matrix33d::Zero();
    H.block<3,3>(3,3) = R.transpose();


    F0 = ((H.transpose()).inverse())*F_end;

    //debug
    std::cout<<"F0 "<<F0<<std::endl;

    F0_after<<F0(0),F0(1),F0(2)-G,F0(3)+G*frame_L.p(1),F0(4)-G*frame_L.p(0),F0(5);

    //debug
    std::cout<<"F0_after"<<F0_after<<std::endl;

    //force to speed and position
    mapping(F0_after);
    now_position<<frame_ur.p(0),frame_ur.p(1),frame_ur.p(2),gammar_ur,beta_ur,alfa_ur;
    des_position = now_position + delta_position;

    //debug
    std::cout<<"delta_position"<<delta_position<<std::endl;
    std::cout<<"des_speed"<<des_speed<<std::endl;

    jointpositions = IKbasedonQP::IK(jointpositions,jointpositions_now,des_position,des_speed,
                           my_chain,bound_low,bound_upper);


    for(unsigned int i=0;i<6;i++){
      jointpositions_msg.position.push_back(jointpositions(i));
    }

    //to real robot
    chatter_pub[6].publish(jointpositions_msg);
    return jointpositions;
//    //debug
//    std_msgs::Float64 msg_pos;
//    msg_pos = des_position[0];
//    chatter_pub[8].publish(msg_pos);


}

/**
 * \brief this function is used for input the waypoint
 * @param double x y z - position of waypoint
 * @param double gammar beta alfa - attitude of waypoint
 * @param t - time stamp
 * @param cartpos_interplotation - used to store the pose of waypoints
 * @param T - used to store the timepoint of waypoints
   */
void input(double x,double y,double z,double gammar, double beta,double alfa,double t,
           std::vector<std::vector<double>> &cartpos_interplotation, std::vector<double> &T){

    std::vector<double> cartpos_pass(6,0);

    cartpos_pass[0]  = x;
    cartpos_pass[1]  = y;
    cartpos_pass[2]  = z;
    cartpos_pass[3]  = gammar;
    cartpos_pass[4]  = beta;
    cartpos_pass[5]  = alfa;

    cartpos_interplotation.push_back(cartpos_pass);
    T.push_back(t);

}


int main(int argc, char **argv)
{

    remove("simul_cartpos.txt");
    remove("cartpos.txt");
    remove("distance.txt");
    remove("cartpos_temp.txt");
    remove("initial eulerangle1.txt");
    remove("des_a.txt");
    remove("initial eulerangle_loop.txt");
    remove("initial eulerangle.txt");
    remove("cartpos_tosend.txt");
    remove("des_cartpos.txt");
    remove("real_cartpos.txt");
    remove("jointpositions after cubic interpolation.txt");
    remove("des_cartendvel.txt");
    ros::init(argc, argv, "torque_teaching");
    ros::NodeHandle n;

    //Publisher
    ros::Publisher chatter_pub1 = n.advertise<std_msgs::Float64>("/shoulder_pan_joint_position_controller/command", 1000);
    ros::Publisher chatter_pub2 = n.advertise<std_msgs::Float64>("/shoulder_lift_joint_position_controller/command", 1000);
    ros::Publisher chatter_pub3 = n.advertise<std_msgs::Float64>("/elbow_joint_position_controller/command", 1000);
    ros::Publisher chatter_pub4 = n.advertise<std_msgs::Float64>("/wrist_1_joint_position_controller/command", 1000);
    ros::Publisher chatter_pub5 = n.advertise<std_msgs::Float64>("/wrist_2_joint_position_controller/command", 1000);
    ros::Publisher chatter_pub6 = n.advertise<std_msgs::Float64>("/wrist_3_joint_position_controller/command", 1000);
    ros::Publisher chatter_pub7 = n.advertise<sensor_msgs::JointState>("/ur5_jointstates",1000);

    //debug
    ros::Publisher des_position_pub = n.advertise<std_msgs::Float64>("/des_position",1000);

    // push back the pub
    std::vector<ros::Publisher> chatter_pub;
    chatter_pub.push_back(chatter_pub1);
    chatter_pub.push_back(chatter_pub2);
    chatter_pub.push_back(chatter_pub3);
    chatter_pub.push_back(chatter_pub4);
    chatter_pub.push_back(chatter_pub5);
    chatter_pub.push_back(chatter_pub6);
    chatter_pub.push_back(chatter_pub7);
    chatter_pub.push_back(des_position_pub);


//Subscriber
    ros::Subscriber sub_realrobot_joint_states = n.subscribe("/joint_states_real_robot",1000,realrobotcallback_JointStates);
    ros::Subscriber sub_gazeborobot_joint_states = n.subscribe("/joint_states",1000,gazeborobotcallback_JointStates);
    ros::Subscriber torque_sensor = n.subscribe("robotiq_force_torque_sensor",100,reCallback);

//time loop
    ros::Rate loop_rate(1/servoj_time);

//jiexi kdl tree
    KDL::Tree my_tree;
    ros::NodeHandle node;
    std::string robot_desc_string;
    node.param("robot_description", robot_desc_string, std::string());
    if (!kdl_parser::treeFromString(robot_desc_string, my_tree)){
        ROS_ERROR("Failed to construct kdl tree");
        return false;
    }

//kdl tree to kdl chain
    KDL::Chain my_chain;
    const std::string &eef = findEndeffector(my_tree);
        if (!my_tree.getChain(my_tree.getRootSegment()->first, eef, my_chain)) {
           ROS_ERROR_STREAM("Could not find chain to " << eef);
           return EXIT_FAILURE;
        }

//    chain_robotiq = chain_ur;
//    KDL::Segment segment_robotiq = chain_robotiq.getSegment(7);

//    segment_robotiq = Segment(Joint(Joint::RotZ),Frame(Rotation::RPY(),Vector()));

//define the passing points
    double alfa;
    double beta;
    double gammar;
    KDL::Frame cartpos;
    KDL::JntArray jointpositions_standard(6);
    jointpositions_standard.data <<-(91.71/180)*PI,(-98.96/180)*PI,(-126.22/180)*PI,(-46.29/180)*PI,(91.39/180)*PI,(-1.78/180)*PI;
    cartpos = IKbasedonQP::FK(jointpositions_standard,my_chain);
    cartpos.M.GetEulerZYX(alfa,beta,gammar);

// calibration points
//1
    std::vector<std::vector<double> > cartpos_interplotation1;
    std::vector<double> T1;
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta,alfa,0,cartpos_interplotation1,T1);
    input(cartpos.p(0)+0.1,cartpos.p(1),cartpos.p(2),gammar+0.3,beta,alfa,1,cartpos_interplotation1,T1);

//2
    std::vector<std::vector<double> > cartpos_interplotation2;
    std::vector<double> T2;
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar+0.3,beta,alfa,2,cartpos_interplotation2,T2);
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar+0.6,beta,alfa,3,cartpos_interplotation2,T2);

//3
    std::vector<std::vector<double> > cartpos_interplotation3;
    std::vector<double> T3;
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar+0.6,beta,alfa,4,cartpos_interplotation3,T3);
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar+0.9,beta,alfa,5,cartpos_interplotation3,T3);

//4
    std::vector<std::vector<double> > cartpos_interplotation4;
    std::vector<double> T4;
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar+0.9,beta,alfa,6,cartpos_interplotation4,T4);
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta+0.3,alfa,7,cartpos_interplotation4,T4);

//5
    std::vector<std::vector<double> > cartpos_interplotation5;
    std::vector<double> T5;
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta+0.3,alfa,8,cartpos_interplotation5,T5);
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta+0.6,alfa,9,cartpos_interplotation5,T5);

//6
    std::vector<std::vector<double> > cartpos_interplotation6;
    std::vector<double> T6;
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta+0.6,alfa,10,cartpos_interplotation6,T6);
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta+0.9,alfa,11,cartpos_interplotation6,T6);

//7
    std::vector<std::vector<double> > cartpos_interplotation7;
    std::vector<double> T7;
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta+0.9,alfa,12,cartpos_interplotation7,T7);
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta+0.2,alfa,13,cartpos_interplotation7,T7);

//8
    std::vector<std::vector<double> > cartpos_interplotation8;
    std::vector<double> T8;
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta+0.2,alfa,14,cartpos_interplotation8,T8);
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta-0.2,alfa,15,cartpos_interplotation8,T8);

//9
    std::vector<std::vector<double> > cartpos_interplotation9;
    std::vector<double> T9;
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta-0.4,alfa,16,cartpos_interplotation9,T9);
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta-0.6,alfa,17,cartpos_interplotation9,T9);

    while (ros::ok()) {

        if (isGetRealJointValue){

            ros::spinOnce();
            jointpositions = homing(jointpositions_now_first,t,0,5,chatter_pub);
            jointpositions = trajectory_do(t,T1,cartpos_interplotation1,my_chain,chatter_pub);
            if(t>T1.front()&&t<(T1.front()+0.003)){
                force_torque_calibration[0] = force_torque;
            }

            jointpositions = trajectory_do(t,T2,cartpos_interplotation2,my_chain,chatter_pub);
            if(t>T2.front()&&t<(T2.front()+0.003)){
                force_torque_calibration[1] = force_torque;

            }

            jointpositions = trajectory_do(t,T3,cartpos_interplotation3,my_chain,chatter_pub);
            if(t>T3.front()&&t<(T3.front()+0.003)){
                force_torque_calibration[2] = force_torque;

            }


            jointpositions = trajectory_do(t,T4,cartpos_interplotation4,my_chain,chatter_pub);
            if(t>T4.front()&&t<(T4.front()+0.003)){
                force_torque_calibration[3] = force_torque;

            }


            jointpositions = trajectory_do(t,T5,cartpos_interplotation5,my_chain,chatter_pub);
            if(t>T5.front()&&t<(T5.front()+0.003)){
                force_torque_calibration[4] = force_torque;

            }

            jointpositions = trajectory_do(t,T6,cartpos_interplotation6,my_chain,chatter_pub);
            if(t>T6.front()&&t<(T6.front()+0.003)){
                force_torque_calibration[5] =force_torque;

            }


            jointpositions = trajectory_do(t,T7,cartpos_interplotation7,my_chain,chatter_pub);
            if(t>T7.front()&&t<(T7.front()+0.003)){
                force_torque_calibration[6] = force_torque;

            }

            jointpositions = trajectory_do(t,T8,cartpos_interplotation8,my_chain,chatter_pub);
            if(t>T8.front()&&t<(T8.front()+0.003)){
                force_torque_calibration[7] = force_torque;

            }


            jointpositions = trajectory_do(t,T9,cartpos_interplotation9,my_chain,chatter_pub);
            if(t>T9.front()&&t<(T9.front()+0.003)){
                force_torque_calibration[8] = force_torque;

            }


            if(t>22&&t<23){
                //debug
                std::cout<<"force_torque_calibration: "<<std::endl;
                std::cout<<force_torque_calibration[0][0]<<std::endl;
                std::cout<<force_torque_calibration[1][0]<<std::endl;
                std::cout<<force_torque_calibration[2][0]<<std::endl;
                std::cout<<force_torque_calibration[3][0]<<std::endl;
                std::cout<<force_torque_calibration[4][0]<<std::endl;
                std::cout<<force_torque_calibration[5][0]<<std::endl;
                std::cout<<force_torque_calibration[6][0]<<std::endl;
                std::cout<<force_torque_calibration[7][0]<<std::endl;
                std::cout<<force_torque_calibration[8][0]<<std::endl;

                centroid_weight_calibration(force_torque_calibration);
            }
            if(t>24){
                force_chasing(jointpositions,my_chain,chatter_pub);
            }

        }//if (isGetRealJointValue)

        else{
            std::cout<<"not get the current jointstates"<<std::endl;
        }// else

        std::cout<<"time:"<<t<<std::endl;
        ros::spinOnce();
        loop_rate.sleep();
        t=t+servoj_time;

    }// while (ros::ok())
    return 0;
}
  
