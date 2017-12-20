#include "ikbasedonqp.h"
#include <Eigen/Dense>
#include "optimal_distance.h"


#define PI 3.14159

USING_NAMESPACE_QPOASES
using namespace std;

extern double servoj_time;
extern std::vector<KDL::Chain> chain_series;
KDL::JntArray joint_in(6);
KDL::JntArray joint_pre(6);
KDL::Jacobian ur_jac(6);
KDL::Jacobian ur_jac_pre(6);

typedef Eigen::Matrix<double,6,6> Matrix6d;
typedef Eigen::Matrix<double,6,5> Matrix65d;
typedef Eigen::Matrix<double,6,1> Matrix61d;
typedef Eigen::Matrix<double,1,6> Matrix16d;
typedef Eigen::Matrix<double,2,6> Matrix26d;
typedef Eigen::Matrix<double,6,Eigen::Dynamic> Matrix6nd;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Matrixxd;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> Matrixx1d;
typedef Eigen::Matrix<double,3,3> Matrix33d;
typedef Eigen::Matrix<double,3,1> Matrix31d;
typedef Eigen::Matrix<double,1,3> Matrix13d;
typedef Eigen::Matrix<double,3,6> Matrix36d;
typedef Eigen::Matrix<double,12,6> Matrix126d;
typedef Eigen::Matrix<double,12,1> Matrix121d;

/**
 * \brief this function is used to get a V shape curve
 * @param t - time stamp
 * @param t_s - start time
 * @param t_f- end time
 * @return the value of cubic interpotation function
 */
double ploy3_we(double t,double t_s,double t_f)
{
    double t_mid=(t_f-t_s)/2;
    double we;
    if(t>=t_s&&t<=t_mid){

      we = 1+(3/((t_mid-t_s)*(t_mid-t_s)))*(0-1)*(t-t_s)*(t-t_s)
          -(2/((t_mid-t_s)*(t_mid-t_s)*(t_mid-t_s)))*(0-1)*(t-t_s)*(t-t_s)*(t-t_s);

      }

    if(t>=t_mid&&t<=t_f){

      we = 0+(3/((t_f-t_mid)*(t_f-t_mid)))*(1-0)*(t-t_mid)*(t-t_mid)
          -(2/((t_f-t_mid)*(t_f-t_mid)*(t_f-t_mid)))*(1-0)*(t-t_mid)*(t-t_mid)*(t-t_mid);

      }

    return we;

}

IKbasedonQP::IKbasedonQP()
{

}

/**
 * \brief this function is  forward kinematics
 * @param jointpositions - joint positions of ur
 * @param my_chain - chain of ur
 * @return the cartpos of ur
 */
KDL::Frame IKbasedonQP::FK(const KDL::JntArray jointpositions,
                           KDL::Chain& my_chain){

    KDL::ChainFkSolverPos_recursive fksolver = KDL::ChainFkSolverPos_recursive(my_chain);
    KDL::Frame cartpos;
    bool kinematics_status;
    kinematics_status = fksolver.JntToCart(jointpositions,cartpos);
        if(kinematics_status>=0){
          printf("%s \n","Succes, thanks KDL!");
          return cartpos;

        }else{
                printf("%s \n","Error: could not calculate forward kinematics :(");
                return cartpos;

        }

}


/**
 * \brief this function is inverse kinematics(no obstacle avoidance)
 * @param jointpositions - joint positions of ur
 * @param my_chain - chain of ur
 * @param jointpositions
 * @param jointpositions_now - current joint position of ur
 * @param des_cartpos - desire cartesian position
 * @param des_cartendvel - desiare cartesian velocity
 * @param bound_low - the low limit of joint
 * @param bound_upper - the upper limit of joint
 * @return the joint position via the IK
 */
KDL::JntArray IKbasedonQP::IK(KDL::JntArray jointpositions,KDL::JntArray jointpositions_now,Matrix61d des_cartpos,
                              Matrix61d des_cartendvel, KDL::Chain& my_chain,
                              Matrixx1d bound_low, Matrixx1d bound_upper){

    Matrix6nd j;
    Matrix61d output_cartendvel;
    Matrix61d jointvel;


    // caculate the Jacobian
    KDL::ChainJntToJacSolver ur_JacSlover = KDL::ChainJntToJacSolver(my_chain);
    unsigned int nj = my_chain.getNrOfJoints();
    KDL::JntArray joint_in(nj);
    KDL::Jacobian ur_jac(nj);

//    joint_in = jointpositions;
    joint_in = jointpositions_now;
    ur_JacSlover.JntToJac(joint_in, ur_jac, -1);
    std::cout<<"Jacobian "<< std::endl << ur_jac.data << std::endl;

  // fk
    Matrix61d real_cartpos;
    KDL::Frame cartpos;
    KDL::Frame simul_cartpos;
    cartpos = IKbasedonQP::FK(jointpositions_now,my_chain);
    simul_cartpos = IKbasedonQP::FK(jointpositions,my_chain);
  //  caculate the Euler

     double alfa;
     double beta;
     double gammar;
     cartpos.M.GetEulerZYX(alfa,beta,gammar);
     real_cartpos << cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta,alfa;


   // estabish the QPmodel

    Matrixxd H1;
    Matrixx1d g1;
    Matrix6d K;
    j = ur_jac.data;
    H1 = j.transpose()*j;
    K << 0.7,0,0,0,0,0,
         0,0.7,0,0,0,0,
         0,0,0.7,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0,
         0,0,0,0,0,0;

    g1 = -j.transpose()*(des_cartendvel+K*(des_cartpos- real_cartpos));
//    g1 = -j.transpose()*(des_cartendvel);
    std::ofstream ofile;


    ofile.open("simul_cartpos.txt",std::ios_base::app);
    ofile<<simul_cartpos.p(0)<<" "
        << simul_cartpos.p(1)<<" "
        << simul_cartpos.p(2)<<" "
        <<std::endl;
    ofile.close();


   //output the real cartendpos

    ofile.open("real_cartpos.txt",std::ios_base::app);
    ofile<< real_cartpos(0)<<" "
        << real_cartpos(1)<<" "
        << real_cartpos(2)<<" "
        << real_cartpos(3)<<" "
        << real_cartpos(4)<<" "
        << real_cartpos(5)<<std::endl;
    ofile.close();

    //output the cartendvel
    ofile.open("des_cartendvel.txt",std::ios_base::app);
    ofile<< des_cartendvel(0)<< " "
        << des_cartendvel(1)<<" "
        << des_cartendvel(2)<<" "
        << des_cartendvel(3)<<" "
        << des_cartendvel(4)<<" "
        << des_cartendvel(5)<<std::endl;
    ofile.close();

    /* Setup data of first QP. */
    real_t H[nj*nj];
    real_t g[nj];

    for(unsigned int i=0;i<nj;i++){
        for(unsigned int j=0;j<nj ;j++){
            H[i*nj+j]=H1(i,j);
        }
    }


    for(unsigned int i=0;i<nj;i++){
        g[i]=g1(i);
    }


    real_t lb[nj];
    real_t ub[nj];

    for(unsigned int i=0;i<nj;i++){
        lb[i] = bound_low(i);
    }
    for(unsigned int i=0;i<nj;i++){
        ub[i] = bound_upper(i);
    }

    /* Setting up QProblemB object. */
    QProblemB example( nj );
    Options options;


    //options.enableFlippingBounds = BT_FALSE;
    options.initialStatusBounds = ST_INACTIVE;
    options.numRefinementSteps = 1;
    options.enableCholeskyRefactorisation = 1;
    example.setOptions( options );

    /* Solve first QP. */
    int_t nWSR = 10;
    example.init( H,g,lb,ub, nWSR,0 );

    /* Get and print solution of first QP. */
    real_t xOpt[nj];
    example.getPrimalSolution( xOpt );


       for(unsigned int i=0;i<nj ;i++){
           jointvel(i)= xOpt[i];
       }


    output_cartendvel = j*jointvel;
    jointpositions.data = jointpositions.data + jointvel*servoj_time;
    return jointpositions;


}

/**
 * \brief this function is  forward kinematics for diffirent chain
 * @param jointpositions - joint positions of ur
 * @return the series cartpos of ur
 */
std::vector<KDL::Frame>  fk_series(KDL::JntArray jointpositions){


    KDL::JntArray jointpositions5(5);
    KDL::JntArray jointpositions4(4);
    KDL::JntArray jointpositions3(3);
    KDL::JntArray jointpositions2(2);

    jointpositions5(0) = jointpositions(0);
    jointpositions5(1) = jointpositions(1);
    jointpositions5(2) = jointpositions(2);
    jointpositions5(3) = jointpositions(3);
    jointpositions5(4) = jointpositions(4);


    std::vector<KDL::Frame> cartpos_series;
    KDL::Frame  cartpos7;
    KDL::Frame  cartpos5;
    KDL::Frame  cartpos4;
    KDL::Frame  cartpos3;
    KDL::Frame  cartpos2;
    cartpos7 = IKbasedonQP::FK(jointpositions,chain_series[6]);
    cartpos5 = IKbasedonQP::FK(jointpositions5,chain_series[5]);
    cartpos4 = IKbasedonQP::FK(jointpositions4,chain_series[4]);
    cartpos3 = IKbasedonQP::FK(jointpositions3,chain_series[3]);
    cartpos2 = IKbasedonQP::FK(jointpositions2,chain_series[2]);
    cartpos_series.push_back(cartpos7);
    cartpos_series.push_back(cartpos5);
    cartpos_series.push_back(cartpos4);
    cartpos_series.push_back(cartpos3);
    cartpos_series.push_back(cartpos2);

    return cartpos_series;
}

/**
 * \brief this function can caculate the jocobi at point of the rigid body
 * @param ur_jac - joint positions of ur
 * @param cartpos - chain of ur
 * @param p - the pose of the point
 * @return the jocobi of the point
 */
Matrix6d jac_atpointp(KDL::Jacobian ur_jac, KDL::Frame cartpos,Matrix31d p){

  if(ur_jac.data.cols()==6){
    std::cout<<"dwadwaijahicahic"<<111111<<std::endl;
      Matrix6d j;
      j = ur_jac.data;
      Matrix31d p0;
      p0<<cartpos.p(0),cartpos.p(1),cartpos.p(2);

      Matrix33d R;
      R<<cartpos.M.data[0],cartpos.M.data[1],cartpos.M.data[2],
          cartpos.M.data[3],cartpos.M.data[4],cartpos.M.data[5],
          cartpos.M.data[6],cartpos.M.data[7],cartpos.M.data[8];
      //get the fanduichengzhen
      Matrix31d z1,z2,z3,z4,z5,z6;
      z1<<j(0,0),j(1,0),j(2,0);
      z2<<j(0,1),j(1,1),j(2,1);
      z3<<j(0,2),j(1,2),j(2,2);
      z4<<j(0,3),j(1,3),j(2,3);
      z5<<j(0,4),j(1,4),j(2,4);
      z6<<j(0,5),j(1,5),j(2,5);

      Matrix33d pz1,pz2,pz3,pz4,pz5,pz6;

      pz1<<0,-z1(2),z1(1),
           z1(2),0,-z1(0),
          -z1(1),z1(0),0;

      pz2<<0,-z2(2),z2(1),
           z2(2),0,-z2(0),
          -z2(1),z2(0),0;

      pz3<<0,-z3(2),z3(1),
           z3(2),0,-z3(0),
          -z3(1),z3(0),0;

      pz4<<0,-z4(2),z4(1),
           z4(2),0,-z4(0),
          -z4(1),z4(0),0;

      pz5<<0,-z5(2),z5(1),
           z5(2),0,-z5(0),
          -z5(1),z5(0),0;

      pz6<<0,-z6(2),z6(1),
           z6(2),0,-z6(0),
          -z6(1),z6(0),0;

    //caculate the jac of the point
     Matrix31d ba,rab;
     ba = p-p0;
     rab= R*ba;

     Matrix31d detaj1,detaj2,detaj3,detaj4,detaj5,detaj6;
     detaj1= pz1*rab;
     detaj2= pz2*rab;
     detaj3= pz3*rab;
     detaj4= pz4*rab;
     detaj5= pz5*rab;
     detaj6= pz6*rab;

     Matrix31d j_new1,j_new2,j_new3,j_new4,j_new5,j_new6;
     j_new1 = z1+detaj1;
     j_new2 = z2+detaj2;
     j_new3 = z3+detaj3;
     j_new4 = z4+detaj4;
     j_new5 = z5+detaj5;
     j_new6 = z6+detaj6;

     Matrix6d j_new;
     j_new<<j_new1(0),j_new2(0),j_new3(0),j_new4(0),j_new5(0),j_new6(0),
            j_new1(1),j_new2(1),j_new3(1),j_new4(1),j_new5(1),j_new6(1),
            j_new1(2),j_new2(2),j_new3(2),j_new4(2),j_new5(2),j_new6(2),
            j(3,0),j(3,1),j(3,2),j(3,3),j(3,4),j(3,5),
            j(4,0),j(4,1),j(4,2),j(4,3),j(4,4),j(4,5),
            j(5,0),j(5,1),j(5,2),j(5,3),j(5,4),j(5,5);
     return j_new;

    //      std::cout << "j_new"<< std::endl<<j_new<<std::endl;
 }
  if(ur_jac.data.cols()==5){
    std::cout<<"dwadwaijahicahic"<<2222222<<std::endl;

    Matrix65d j;
    j = ur_jac.data;
    Matrix31d p0;
    p0<<cartpos.p(0),cartpos.p(1),cartpos.p(2);

    Matrix33d R;
    R<<cartpos.M.data[0],cartpos.M.data[1],cartpos.M.data[2],
        cartpos.M.data[3],cartpos.M.data[4],cartpos.M.data[5],
        cartpos.M.data[6],cartpos.M.data[7],cartpos.M.data[8];
    //get the fanduichengzhen
    Matrix31d z1,z2,z3,z4,z5;
    z1<<j(0,0),j(1,0),j(2,0);
    z2<<j(0,1),j(1,1),j(2,1);
    z3<<j(0,2),j(1,2),j(2,2);
    z4<<j(0,3),j(1,3),j(2,3);
    z5<<j(0,4),j(1,4),j(2,4);

    Matrix33d pz1,pz2,pz3,pz4,pz5;

    pz1<<0,-z1(2),z1(1),
         z1(2),0,-z1(0),
        -z1(1),z1(0),0;

    pz2<<0,-z2(2),z2(1),
         z2(2),0,-z2(0),
        -z2(1),z2(0),0;

    pz3<<0,-z3(2),z3(1),
         z3(2),0,-z3(0),
        -z3(1),z3(0),0;

    pz4<<0,-z4(2),z4(1),
         z4(2),0,-z4(0),
        -z4(1),z4(0),0;

    pz5<<0,-z5(2),z5(1),
         z5(2),0,-z5(0),
        -z5(1),z5(0),0;

   //caculate the jac of the point
   Matrix31d ba,rab;
   ba = p-p0;
   rab= R*ba;

   Matrix31d detaj1,detaj2,detaj3,detaj4,detaj5;
   detaj1= pz1*rab;
   detaj2= pz2*rab;
   detaj3= pz3*rab;
   detaj4= pz4*rab;
   detaj5= pz5*rab;

   Matrix31d j_new1,j_new2,j_new3,j_new4,j_new5;
   j_new1 = z1+detaj1;
   j_new2 = z2+detaj2;
   j_new3 = z3+detaj3;
   j_new4 = z4+detaj4;
   j_new5 = z5+detaj5;

   Matrix6d j_new;
   j_new<<j_new1(0),j_new2(0),j_new3(0),j_new4(0),j_new5(0),0,
          j_new1(1),j_new2(1),j_new3(1),j_new4(1),j_new5(1),0,
          j_new1(2),j_new2(2),j_new3(2),j_new4(2),j_new5(2),0,
          j(3,0),j(3,1),j(3,2),j(3,3),j(3,4),0,
          j(4,0),j(4,1),j(4,2),j(4,3),j(4,4),0,
          j(5,0),j(5,1),j(5,2),j(5,3),j(5,4),0;
   return j_new;
  }


}

/**
 * \brief this function is inverse kinematics(obstacle avoidance)
 * @param jointpositions - joint positions of ur
 * @param my_chain - chain of ur
 * @param jointpositions
 * @param jointpositions_now - current joint position of ur
 * @param des_cartpos - desire cartesian position
 * @param des_cartendvel - desiare cartesian velocity
 * @param bound_low - the low limit of joint
 * @param bound_upper - the upper limit of joint
 * @param t - time stampe
 * @param T - used to store the timepoint
 * @return the joint position via the IK
 */

KDL::JntArray IKbasedonQP::IK(KDL::JntArray jointpositions,KDL::JntArray jointpositions_now,Matrix61d des_cartpos,
                              Matrix61d des_cartendvel, KDL::Chain& my_chain,
                              Matrixx1d bound_low, Matrixx1d bound_upper,std::vector<double> T,double t){

    std::ofstream ofile;
    Matrix61d output_cartendvel;
    Matrix61d jointvel;

    // caculate the Jacobian end
    KDL::ChainJntToJacSolver ur_JacSlover = KDL::ChainJntToJacSolver(my_chain);
    unsigned int nj = my_chain.getNrOfJoints();



    joint_in = jointpositions_now;


        joint_pre(0) = joint_in(0)-0.01;
        joint_pre(1) = joint_in(1)-0.01;
        joint_pre(2) = joint_in(2)-0.01;
        joint_pre(3) = joint_in(3)-0.01;
        joint_pre(4) = joint_in(4)-0.01;
        joint_pre(5) = joint_in(5)-0.01;




    ur_JacSlover.JntToJac(joint_in, ur_jac, -1);
    ur_JacSlover.JntToJac(joint_pre, ur_jac_pre, -1);


    std::cout<<"Jacobian "<< std::endl << ur_jac.data << std::endl;
    std::cout<<"Jacobian pre "<< std::endl << ur_jac_pre.data << std::endl;
    Matrix6d j = ur_jac.data;
    Matrix6d j_pre = ur_jac_pre.data;


    // caculate the Jacobian joint5
    KDL::ChainJntToJacSolver ur_JacSlover5 = KDL::ChainJntToJacSolver(chain_series[5]);
    KDL::JntArray joint5(5);
    KDL::Jacobian ur_jac5(5);
    joint5(0) = jointpositions_now(0);
    joint5(1) = jointpositions_now(1);
    joint5(2) = jointpositions_now(2);
    joint5(3) = jointpositions_now(3);
    joint5(4) = jointpositions_now(4);
    ur_JacSlover5.JntToJac(joint5, ur_jac5, -1);
    std::cout<<"Jacobian5 "<< std::endl << ur_jac5.data << std::endl;


   // fk
    Matrix61d real_cartpos;
    std::vector<KDL::Frame> cartpos_series;
    KDL::Frame  cartpos;
    cartpos_series = fk_series(jointpositions_now);
    cartpos = cartpos_series[0];
    double alfa;
    double beta;
    double gammar;
    cartpos.M.GetEulerZYX(alfa,beta,gammar);
    std::cout<<"cartpos"<<std::endl<<cartpos<<std::endl;


    KDL::Frame simul_cartpos;
    simul_cartpos = FK(jointpositions,my_chain);
    ofile.open("simul_cartpos.txt",std::ios_base::app);
    ofile<<simul_cartpos.p(0)<<" "
        << simul_cartpos.p(1)<<" "
        << simul_cartpos.p(2)<<" "
        <<std::endl;
    ofile.close();

    //fk5
    KDL::Frame  cartpos5;
    cartpos5 = cartpos_series[1];
    double alfa5;
    double beta5;
    double gammar5;
    cartpos5.M.GetEulerZYX(alfa5,beta5,gammar5);
    std::cout<<"cartpos5"<<std::endl<<cartpos5<<std::endl;

    //obstacle location fk
    double alfa_obs;
    double beta_obs;
    double gammar_obs;
    KDL::JntArray jointpositions_standard(6);
    jointpositions_standard.data <<-(91.71/180)*PI,(-98.96/180)*PI,(-126.22/180)*PI,(-46.29/180)*PI,(91.39/180)*PI,(-1.78/180)*PI;
    KDL::Frame cartpos_obs_fk;
    cartpos_obs_fk = IKbasedonQP::FK(jointpositions_standard,my_chain);
    cartpos_obs_fk.M.GetEulerZYX(alfa_obs,beta_obs,gammar_obs);

    //pair1 end and obstacle
     Matrix31d p;
     Matrix31d q;
     std::vector<Matrix31d> geometric_info(2);

    //cacaulate the realtime distance
     double d;
     geometric_info= od(cartpos,d);
     std::cout<<"d"<<d<<std::endl;
     p = geometric_info[0];
     q = geometric_info[1];

     Matrix61d pq6,p6,q6;
     p6<<p(0),p(1),p(2),gammar,beta,alfa;
     q6<<q(0),q(1),q(2),gammar_obs,beta_obs,alfa_obs;
     pq6 = p6-q6;
     double d1;
     d1= pq6.norm();
  
     double e=0.5;
     double ds=0.05;
     double di=0.08;
     Matrix16d A1;
     Matrix16d nn;
     nn = pq6/d1;
     std::cout<<"nn"<<nn<<std::endl;

     Matrix6d j_new = jac_atpointp(ur_jac,cartpos,p);
     A1= nn*j_new;
     std::cout<<"A1"<<A1<<std::endl;


     //pair2 joint5 and obstacle
     double d5;
     Matrix31d p5;
     Matrix31d q5;
     std::vector<Matrix31d> geometric_info1(2);
     geometric_info1= od1(cartpos5,d5);
     std::cout<<"d5"<<d5<<std::endl;
     p5 = geometric_info1[0];
     q5 = geometric_info1[1];

     Matrix61d pq65,p65,q65;
     p65<<p5(0),p5(1),p5(2),gammar5,beta5,alfa5;
     q65<<q5(0),q5(1),q5(2),gammar_obs,beta_obs,alfa_obs;
     pq65 = p65-q65;
     double len5;
     len5 = pq65.norm();

     double e5=0.5;
     double ds5=0.05;
     double di5=0.1;
     Matrix16d n5;
     n5 = pq65/len5;

     std::cout<<"n5"<<n5<<std::endl;

     Matrix6d j_new5;
     j_new5 = jac_atpointp(ur_jac5,cartpos5,p5);
     Matrix16d A2;
     A2 = n5*j_new5;
     std::cout<<"A2"<<A2<<std::endl;

    Matrix26d AA;
    AA<<A1(0),A1(1),A1(2),A1(3),A1(4),A1(5),
        A2(0),A2(1),A2(2),A2(3),A2(4),A2(5);
    std::cout<<"AA"<<AA<<std::endl;

    // estabish the QPmodel

    //joint tasks
    Matrix6d A_joint;
    Matrix61d b_joint;
    A_joint<< 0.005,0,0,0,0,0,
              0,0.005,0,0,0,0,
              0,0,0.005,0,0,0,
              0,0,0,0.005,0,0,
              0,0,0,0,0.005,0,
              0,0,0,0,0,0.005;
    double mani_now;
    double mani_pre;
    mani_now = sqrt((j*j.transpose()).determinant());
    std::cout<<"mani_now"<<mani_now<<std::endl;

    mani_pre = sqrt((j_pre*j_pre.transpose()).determinant());
    std::cout<<"mani_pre"<<mani_pre<<std::endl;
    b_joint<< (mani_now-mani_pre),
              (mani_now-mani_pre),
              (mani_now-mani_pre),
              (mani_now-mani_pre),
              (mani_now-mani_pre),
              (mani_now-mani_pre);

    //cartesian tasks

    double we;
    we = ploy3_we(t,T.front(),T.back());
    Matrix6d We;
    We << we,0,0,0,0,0,
          0,we,0,0,0,0,
          0,0,we,0,0,0,
          0,0,0,we,0,0,
          0,0,0,0,we,0,
          0,0,0,0,0,we;


    Matrix6d H1;
    Matrix61d g1;
    Matrix6d K;
    Matrix6d A_cart;
    Matrix126d A_total;
    Matrix121d b_total;
    A_cart =We*j;
    K << 0.7,0,0,0,0,0,
        0,0.7,0,0,0,0,
        0,0,0.7,0,0,0,
        0,0,0,0,0,0,
        0,0,0,0,0,0,
        0,0,0,0,0,0;

    real_cartpos << cartpos.p(0),cartpos.p(1),cartpos.p(2),gammar,beta,alfa;
    Matrix61d b_cart;
    b_cart =des_cartendvel+K*(des_cartpos- real_cartpos);
    A_total<<A_cart,A_joint;
    b_total<<b_cart,b_joint;
    std::cout<<"A_total"<<A_total<<std::endl;
    std::cout<<"b_total"<<b_total<<std::endl;
    g1 = -A_total.transpose()*b_total;
    H1 = A_total.transpose()*A_total;

//    g1 = -A_cart.transpose()*b_cart;
//    H1 = A_cart.transpose()*A_cart;

    //    g1 = -j.transpose()*(des_cartendvel);


   //output the desired cartendpos

    ofile.open("des_cartpos.txt",std::ios_base::app);
    ofile<< des_cartpos(0)<<" "
        << des_cartpos(1)<<" "
        << des_cartpos(2)<<" "
        << des_cartpos(3)<<" "
        << des_cartpos(4)<<" "
        << des_cartpos(5)<<std::endl;
    ofile.close();

   //output the real cartendpos

    ofile.open("real_cartpos.txt",std::ios_base::app);
    ofile<< real_cartpos(0)<<" "
        << real_cartpos(1)<<" "
        << real_cartpos(2)<<" "
        << real_cartpos(3)<<" "
        << real_cartpos(4)<<" "
        << real_cartpos(5)<<std::endl;
    ofile.close();

    //output the cartendvel
    ofile.open("des_cartendvel.txt",std::ios_base::app);
    ofile<< des_cartendvel(0)<< " "
        << des_cartendvel(1)<<" "
        << des_cartendvel(2)<<" "
        << des_cartendvel(3)<<" "
        << des_cartendvel(4)<<" "
        << des_cartendvel(5)<<std::endl;
    ofile.close();

    /* Setup data of first QP. */
    real_t H[nj*nj];
    real_t A[2*6];
//    real_t A[6];
    real_t g[nj];
    real_t lbA[2] = { -e*(d-ds)/(di-ds),-e5*(d5-ds5)/(di5-ds5)};
    real_t ubA[2] = { 100000,100000};
    //***********debug***************//
//    real_t lbA[1] = { -e*(d-ds)/(di-ds)};
//    real_t lbA[1] = { -e5*(d5-ds5)/(di5-ds5)};

//    real_t ubA[1] = { 100000};
    //***********debug***************//
    real_t lb[nj];
    real_t ub[nj];

    for(unsigned int i=0;i<nj;i++){
        for(unsigned int j=0;j<nj ;j++){
            H[i*nj+j]=H1(i,j);
        }
    }

    for(unsigned int i=0;i<nj;i++){
        g[i]=g1(i);
    }

    for(unsigned int i=0;i<2;i++){
       for(unsigned int j=0;j<nj ;j++){
         A[i*6+j]=AA(i,j);
      }
    }
  //***********debug***************//
//    for(unsigned int i=0;i<nj;i++){
//          A[i]=A2(i);
//      }
//***********debug***************//

    for(unsigned int i=0;i<nj;i++){
        lb[i] = bound_low(i);
    }
    for(unsigned int i=0;i<nj;i++){
        ub[i] = bound_upper(i);
    }

  /* Setting up SQProblem object. */
    SQProblem example( 6,2 );

    /* Solve first QP. */
    int_t nWSR = 10;
    example.init( H,g,A,lb,ub,lbA,ubA, nWSR,0 );


    /* Get and print solution of first QP. */
    real_t xOpt[6];
    example.getPrimalSolution( xOpt );

    for(unsigned int i=0;i<nj ;i++){
        jointvel(i)= xOpt[i];
    }
    output_cartendvel = j*jointvel;
    jointpositions.data = jointpositions.data + jointvel*servoj_time;
    return jointpositions;

}





