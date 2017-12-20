#ifndef IKBASEDONQP_H
#define IKBASEDONQP_H

#include <Eigen/Dense>
#include <kdl/chain.hpp>
#include <kdl/chainfksolver.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/frames_io.hpp>
#include <kdl/kdl.hpp>
#include <kdl/jntarray.hpp>
#include <qpOASES.hpp>

typedef Eigen::Matrix<double,6,6> Matrix6d;
typedef Eigen::Matrix<double,6,1> Matrix61d;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> Matrixx1d;
typedef Eigen::Matrix<double,4,4> Matrix44d;


class IKbasedonQP
{
public:
  explicit IKbasedonQP();
  ~IKbasedonQP();



  static KDL::Frame FK(const KDL::JntArray jointpositions,KDL::Chain& my_chain);

  static KDL::JntArray IK( KDL::JntArray jointpositions,KDL::JntArray jointpositions_now, Matrix61d des_cartpos,
                           Matrix61d des_cartendvel,KDL::Chain& my_chain,
                           Matrixx1d bound_low,Matrixx1d bound_upper);


};

#endif // IKBASEDONQP_H
