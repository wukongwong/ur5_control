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
#include <vector>
#include <math.h>
#include "ikbasedonqp.h"
#include <SPLINTER/bspline.h>
#include <SPLINTER/datatable.h>
#include <SPLINTER/bsplinebasis.h>
#include "trajectory_generate.h"
#include <kdl/joint.hpp>
using namespace KDL;
using namespace SPLINTER;

#define PI 3.14159

KDL::JntArray jointpositions_now(6);
KDL::JntArray jointpositions_now_first(6);
KDL::JntArray jointpositions_gazebo(6);
std::vector<KDL::Chain> chain_series;
KDL::Frame  cartpos_dh;
KDL::JntArray jointpositions(6);
unsigned int j,m=0;
double alfa_temp;
double beta_temp;
double gammar_temp;
KDL::Frame cartpos_temp;
bool isGetRealJointValue = false;
typedef Eigen::Matrix<double,6,6> Matrix6d;
typedef Eigen::Matrix<double,6,1> Matrix61d;
double servoj_time = 0.003;
double t_homing = 5;
double t=0;



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


/**
 * \brief the cubic interpolation in joint space,move the robot to the standard position from current position
 * @param jointpositions_s - the start jointposition
 * @param t_s - start time
 * @param t_f - end time
 * @param chatter_pub publisher
 * @return  the jointpositions at t time after the cubic interpolation
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


/**
 * \brief generete the planning trajectory and publish the jointposition to realrobot and gazebo
 * @param T_user - time point of user planning
 * @param cartpospoint - waypoint of user planning
 * @param my_chain - chain
 * @param chatter_pub publisher
 * @return  the jointpositions at t time after the trajectory planning
 */
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

//            jointpositions = IKbasedonQP::IK(jointpositions,jointpositions_now,des_cartpos,des_cartendvel,
//                                   my_chain,bound_low,bound_upper,T,t);
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

//            jointpositions = IKbasedonQP::IK(jointpositions,jointpositions_now,des_cartpos,des_cartendvel,
//                                   my_chain,bound_low,bound_upper,T,t);
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

    if(t>=T.back()&&t<=T.back()+2){

        std::cout<<"ok12222"<<std::endl;
        std_msgs::Float64 msg1;
        std_msgs::Float64 msg2;
        std_msgs::Float64 msg3;
        std_msgs::Float64 msg4;
        std_msgs::Float64 msg5;
        std_msgs::Float64 msg6;
        sensor_msgs::JointState jointpositions_msg;

       //two points to three points
        std::vector<double> T_new;
        std::vector<std::vector<double> > cartpospoint_new;
        double T_mid = (T.back()+(T.back()+2))/2;

        T_new.push_back(T.back());
        T_new.push_back(T_mid);
        T_new.push_back(T.back()+2);

        ++m;

        if(m==1){
            cartpos_temp = IKbasedonQP::FK(jointpositions_now,my_chain);
            cartpos_temp.M.GetEulerZYX(alfa_temp,beta_temp,gammar_temp);
         }

//        std::ofstream ofile;
        std::vector<double> pushdata_temp(6,0);
        pushdata_temp[0] = cartpos_temp.p(0);
        pushdata_temp[1] = cartpos_temp.p(1);
        pushdata_temp[2] = cartpos_temp.p(2);
        pushdata_temp[3] = gammar_temp;
        pushdata_temp[4] = beta_temp;
        pushdata_temp[5] = alfa_temp;

        std::cout<<"pushdata_temp[0]"<<pushdata_temp[0]<<std::endl;
        std::vector<std::vector<double> > cartpospoint_temp(2,pushdata_temp);
        cartpospoint_temp[0] = pushdata_temp;
        cartpospoint_temp[1] = cartpospoint.back();

        std::cout<<"cartpospoint[1][0]"<<cartpospoint[1][0]<<std::endl;

        std::vector<double> cartpos_mid(6,0);
        cartpos_mid[0]=(cartpospoint_temp[0][0]+cartpospoint_temp[1][0])/2;
        cartpos_mid[1]=(cartpospoint_temp[0][1]+cartpospoint_temp[1][1])/2;
        cartpos_mid[2]=(cartpospoint_temp[0][2]+cartpospoint_temp[1][2])/2;
        cartpos_mid[3]=(cartpospoint_temp[0][3]+cartpospoint_temp[1][3])/2;
        cartpos_mid[4]=(cartpospoint_temp[0][4]+cartpospoint_temp[1][4])/2;
        cartpos_mid[5]=(cartpospoint_temp[0][5]+cartpospoint_temp[1][5])/2;

        cartpospoint_new.push_back(cartpospoint_temp[0]);
        cartpospoint_new.push_back(cartpos_mid);
        cartpospoint_new.push_back(cartpospoint_temp[1]);
{
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

     //IK

        jointpositions = IKbasedonQP::IK(jointpositions,jointpositions_now,des_cartpos,des_cartendvel,
                               my_chain,bound_low,bound_upper);


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



    }

    return jointpositions;

}


/**
 * \brief make the robot sleep for a certain time
 * @param t_now - start time
 * @param sleeptime - sleep duration
 */
void trajectory_sleep(double t,double t_now,double sleeptime){

    double t_start;
    double t_end;
    t_start= t_now;
    t_end = t_start+sleeptime;
    if(t>=t_start&&t<=t_end){

         ROS_INFO("sleep");

       }
}

// t1 t2 is system time!
/**
 * \brief make the robot do a line
 * @param cartpos_start - the start point of line
 * @param cartpos_end - the end point of line
 * @param t - time stamp
 * @param t1 - start time
 * @param t2 - end time
 * @param my_chain - kdl chain
 * @param chatter_pub - ros publisher
 * @return the jointpositions at t time after the trajectory planning(line)
 */
KDL::JntArray line(std::vector<double> cartpos_start,std::vector<double> cartpos_end,
                   double t, double t1,double t2,int n,KDL::Chain my_chain,
                   std::vector<ros::Publisher> chatter_pub){


    if(t>=t1&&t<=t2){

        std_msgs::Float64 msg1;
        std_msgs::Float64 msg2;
        std_msgs::Float64 msg3;
        std_msgs::Float64 msg4;
        std_msgs::Float64 msg5;
        std_msgs::Float64 msg6;
        sensor_msgs::JointState jointpositions_msg;

        std::vector<double> T;
        std::vector<std::vector<double> > cartpos_inter(n,std::vector<double>(6,0));
//        std::vector<double> deltacartpos(n,0);
        std::vector<double> delta;

        for(unsigned int i=0;i<6;i++){

          delta.push_back((cartpos_end[i]-cartpos_start[i])/(n-1));
        }

         for(unsigned int i=0;i<n;i++){




           cartpos_inter[i][0]=cartpos_start[0]+i*delta[0];
           cartpos_inter[i][1]=cartpos_start[1]+i*delta[1];
           cartpos_inter[i][2]=cartpos_start[2]+i*delta[2];
           cartpos_inter[i][3]=cartpos_start[3]+i*delta[3];
           cartpos_inter[i][4]=cartpos_start[4]+i*delta[4];
           cartpos_inter[i][5]=cartpos_start[5]+i*delta[5];
         }



        double deltat = (t2-t1)/(n-1);

        for(unsigned int i=0;i<n;i++){

          T.push_back(t1+i*deltat);

        }


        trajectory_generate tg;
        tg.spline_generate(T,cartpos_inter);

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
                               my_chain,bound_low,bound_upper,T,t);


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


/**
 * \brief this function is used to get the series of chain.
 * @param my_chain - UR_chain

   */
void chain_get(KDL::Chain my_chain){

    const KDL::Segment &segment0 = my_chain.getSegment(0);
    const KDL::Segment &segment1 = my_chain.getSegment(1);
    const KDL::Segment &segment2 = my_chain.getSegment(2);
    const KDL::Segment &segment3 = my_chain.getSegment(3);
    const KDL::Segment &segment4 = my_chain.getSegment(4);
    const KDL::Segment &segment5 = my_chain.getSegment(5);
    const KDL::Segment &segment6 = my_chain.getSegment(6);
    const KDL::Segment &segment7 = my_chain.getSegment(7);

    KDL::Chain my_chain1;
    my_chain1.addSegment(segment0);

    KDL::Chain my_chain2;
    my_chain2.addSegment(segment0);
    my_chain2.addSegment(segment1);

    KDL::Chain my_chain3;
    my_chain3.addSegment(segment0);
    my_chain3.addSegment(segment1);
    my_chain3.addSegment(segment2);

    KDL::Chain my_chain4;
    my_chain4.addSegment(segment0);
    my_chain4.addSegment(segment1);
    my_chain4.addSegment(segment2);
    my_chain4.addSegment(segment3);

    KDL::Chain my_chain5;
    my_chain5.addSegment(segment0);
    my_chain5.addSegment(segment1);
    my_chain5.addSegment(segment2);
    my_chain5.addSegment(segment3);
    my_chain5.addSegment(segment4);
    my_chain5.addSegment(segment5);

    KDL::Chain my_chain6;
    my_chain6.addSegment(segment0);
    my_chain6.addSegment(segment1);
    my_chain6.addSegment(segment2);
    my_chain6.addSegment(segment3);
    my_chain6.addSegment(segment4);
    my_chain6.addSegment(segment5);
    my_chain6.addSegment(segment6);
    my_chain6.addSegment(segment7);


    chain_series.push_back(my_chain);
    chain_series.push_back(my_chain1);
    chain_series.push_back(my_chain2);
    chain_series.push_back(my_chain3);
    chain_series.push_back(my_chain4);
    chain_series.push_back(my_chain5);
    chain_series.push_back(my_chain6);
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
    ros::init(argc, argv, "talker");
    ros::NodeHandle n;

//Publisher
    ros::Publisher chatter_pub1 = n.advertise<std_msgs::Float64>("/shoulder_pan_joint_position_controller/command", 1000);
    ros::Publisher chatter_pub2 = n.advertise<std_msgs::Float64>("/shoulder_lift_joint_position_controller/command", 1000);
    ros::Publisher chatter_pub3 = n.advertise<std_msgs::Float64>("/elbow_joint_position_controller/command", 1000);
    ros::Publisher chatter_pub4 = n.advertise<std_msgs::Float64>("/wrist_1_joint_position_controller/command", 1000);
    ros::Publisher chatter_pub5 = n.advertise<std_msgs::Float64>("/wrist_2_joint_position_controller/command", 1000);
    ros::Publisher chatter_pub6 = n.advertise<std_msgs::Float64>("/wrist_3_joint_position_controller/command", 1000);
    ros::Publisher chatter_pub7 = n.advertise<sensor_msgs::JointState>("/ur5_jointstates",1000);

// push back the pub
    std::vector<ros::Publisher> chatter_pub;
    chatter_pub.push_back(chatter_pub1);
    chatter_pub.push_back(chatter_pub2);
    chatter_pub.push_back(chatter_pub3);
    chatter_pub.push_back(chatter_pub4);
    chatter_pub.push_back(chatter_pub5);
    chatter_pub.push_back(chatter_pub6);
    chatter_pub.push_back(chatter_pub7);

//Subscriber
    ros::Subscriber sub_realrobot_joint_states = n.subscribe("/joint_states_real_robot",1000,realrobotcallback_JointStates);
    ros::Subscriber sub_gazeborobot_joint_states = n.subscribe("/joint_states",1000,gazeborobotcallback_JointStates);

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
//get the chain_series
   chain_get(my_chain);

// input the passingpoint
    KDL::Frame cartpos;
    double alfa;
    double beta;
    double gammar;
    KDL::JntArray jointpositions_standard(6);
    jointpositions_standard.data <<-(91.71/180)*PI,(-98.96/180)*PI,(-126.22/180)*PI,(-46.29/180)*PI,(91.39/180)*PI,(-1.78/180)*PI;

    cartpos = IKbasedonQP::FK(jointpositions_standard,my_chain);
    cartpos.M.GetEulerZYX(alfa,beta,gammar);

// debug
    std::ofstream ofile;

//output the standard cartpos
    ofile.open("the standard cartpos.txt",std::ios_base::app);
    ofile<<cartpos.p(0)<<" "
         <<cartpos.p(1)<<" "
         <<cartpos.p(2)<<" "
         <<gammar<<" "
         <<beta<<" "
         <<alfa<<" "<<std::endl;

    ofile.close();

// iput the line startpoint and endpoint
    std::vector<std::vector<double>> cartpos_interplotation1;
    std::vector<double> T1;

//point1 and point 2
    input(cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta,alfa,0,cartpos_interplotation1,T1);
    input(cartpos.p(0)+0.3,cartpos.p(1)-0.2,cartpos.p(2),gammar,beta,alfa,2,cartpos_interplotation1,T1);

    while (ros::ok()) {

        if (isGetRealJointValue){

            ros::spinOnce();
            jointpositions = homing(jointpositions_now_first,t,0,5,chatter_pub);
            jointpositions = trajectory_do(t,T1,cartpos_interplotation1,my_chain,chatter_pub);

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
  
