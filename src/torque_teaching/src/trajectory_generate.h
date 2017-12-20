#ifndef TRAJECTORY_GENERATE_H
#define TRAJECTORY_GENERATE_H
#include "spline.h"
#include <vector>
#include <Eigen/Dense>
#include <kdl/chain.hpp>
#include <kdl/chainfksolver.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/frames_io.hpp>
#include <kdl/kdl.hpp>
#include <kdl/jntarray.hpp>

using namespace tk;
class spline;
class trajectory_generate
{
public:

    trajectory_generate();


    std::vector<double> x,y,z,gammar,beta,alfa;

    tk::spline s_x;
    tk::spline s_y;
    tk::spline s_z;
    tk::spline s_gammar;
    tk::spline s_beta;
    tk::spline s_alfa;


    void spline_generate(std::vector<double> t,std::vector<std::vector<double> > cartpos);

    KDL::JntArray homing(KDL::JntArray jointpositions_s,
                         double t,double t_s,double t_f);


};

#endif // TRAJECTORY_GENERATE_H
