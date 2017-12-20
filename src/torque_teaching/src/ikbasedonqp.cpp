#include "ikbasedonqp.h"
#include <Eigen/Dense>


#define PI 3.14159

USING_NAMESPACE_QPOASES
using namespace std;

extern double servoj_time;
KDL::JntArray joint_in(6);
KDL::JntArray joint_pre(6);
KDL::Jacobian ur_jac(6);
KDL::Jacobian ur_jac_pre(6);


typedef Eigen::Matrix<double,4,4> Matrix44d;
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






