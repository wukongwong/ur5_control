#include "serial.h"
#include "ros/ros.h"
#include "std_msgs/Float64.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <angles/angles.h>
#include <Eigen/Dense>
#include <boost/make_shared.hpp>
#include <math.h>
#include <boost/thread.hpp>

unsigned char DataUse1[5];
unsigned char DataUse2[5];



int main(int argc, char **argv)
{

    ros::init(argc, argv, "serial_port");
    ros::NodeHandle n;
//    boost::thread

//Publisher
    ros::Publisher sensordata1 = n.advertise<std_msgs::Float64>("sensordata1", 1000);
    ros::Publisher sensordata2 = n.advertise<std_msgs::Float64>("sensordata2", 1000);
    std_msgs::Float64 sensordata1_pub;
    std_msgs::Float64 sensordata2_pub;

//time loop
    ros::Rate loop_rate(10);

    int tty_fd = tty_open485("/dev/ttyUSB0");
    assert(tty_fd>0);
    // clear zero
//    unsigned char clearzero[6] = {0x01,0x06,0x00,0x60,0x00,0x01};
//    CRC(clearzero,6);
    unsigned char clearzero_crc1[8]={0x01,0x06,0x00,0x60,0x00,0x01,0x48,0x14};
    unsigned char clearzero_crc2[8]={0x02,0x06,0x00,0x60,0x00,0x01,0x48,0x27};
    //send read
    unsigned char read_crc1[8]={0x01,0x03,0x00,0x00,0x00,0x01,0x84,0x0A};
    unsigned char read_crc2[8]={0x02,0x03,0x00,0x00,0x00,0x01,0x84,0x39};

    write(tty_fd,clearzero_crc1,8);
    write(tty_fd,clearzero_crc2,8);
//    float sensordata1;
//    float sensordata2;

    while(ros::ok()){

    //sensor1
        write(tty_fd,read_crc1,8);
        ros::Duration(0.025).sleep(); // sleep for 25ms

        //sensor1
        Data485_Saving(tty_fd,0x01,DataUse1);
        unsigned short data1;
        data1 = (DataUse1[3]<<8)+DataUse1[4];
        printf("data1=%d\n",data1);

        if(data1<6000){
            sensordata1_pub.data = (double)data1/100;
            printf("sensordata1_pub.data=%f\n",sensordata1_pub.data);

            sensordata1.publish(sensordata1_pub);
            std::cout<<"sensordata1 "<<sensordata1_pub.data<<std::endl;
        }
        else{
            sensordata1_pub.data  =(double) (data1-0xFFFF)/100;
             printf("sensordata1_pub.data=%f\n",sensordata1_pub.data);
            sensordata1.publish(sensordata1_pub);
            std::cout<<"sensordata1"<<sensordata1_pub.data<<std::endl;
        }

        //sensor2
        write(tty_fd,read_crc2,8);
        ros::Duration(0.025).sleep(); // sleep for 25ms
        Data485_Saving(tty_fd,0x02,DataUse2);
        unsigned short data2;
        data2 = (DataUse2[3]<<8)+DataUse2[4];
        printf("data2=%d\n",data2);

        if(data2<6000){
            sensordata2_pub.data = (double)data2/100;
            printf("sensordata2_pub.data=%f\n",sensordata2_pub.data);

            sensordata2.publish(sensordata2_pub);
            std::cout<<"sensordata2 "<<sensordata2_pub.data<<std::endl;
        }
        else{
            sensordata2_pub.data  =(double) (data2-0xFFFF)/100;
             printf("sensordata2_pub.data=%f\n",sensordata2_pub.data);
            sensordata2.publish(sensordata2_pub);
            std::cout<<"sensordata2"<<sensordata2_pub.data<<std::endl;

        }
    ros::spinOnce();
    loop_rate.sleep();

    }

    return 0;
}



