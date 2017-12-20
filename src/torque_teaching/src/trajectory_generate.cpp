#include "trajectory_generate.h"

#define PI 3.14159

trajectory_generate::trajectory_generate()
{

}

void trajectory_generate::spline_generate(std::vector<double> t,std::vector<std::vector<double>> cartpos){

    int n= t.size();

    std::vector<double> cartpos_x;
    std::vector<double> cartpos_y;
    std::vector<double> cartpos_z;
    std::vector<double> cartpos_gammar;
    std::vector<double> cartpos_beta;
    std::vector<double> cartpos_alfa;

    // get the interplotate point
        for(unsigned int i=0;i<n;i++){

            cartpos_x.push_back(cartpos[i][0]);

        }
        for(unsigned int i=0;i<n;i++){

            cartpos_y.push_back(cartpos[i][1]) ;

        }
        for(unsigned int i=0;i<n;i++){

            cartpos_z.push_back(cartpos[i][2]); ;

        }
        for(unsigned int i=0;i<n;i++){

            cartpos_gammar.push_back(cartpos[i][3]);

        }
        for(unsigned int i=0;i<n;i++){

            cartpos_beta.push_back(cartpos[i][4]);

        }
        for(unsigned int i=0;i<n;i++){

            cartpos_alfa.push_back(cartpos[i][5]);

        }


    // set the boundary conditions

    s_x.set_boundary(tk::spline::first_deriv,0.0,tk::spline::first_deriv,0,false);
    s_y.set_boundary(tk::spline::first_deriv,0.0,tk::spline::first_deriv,0,false);
    s_z.set_boundary(tk::spline::first_deriv,0.0,tk::spline::first_deriv,0,false);
    s_gammar.set_boundary(tk::spline::first_deriv,0.0,tk::spline::first_deriv,0,false);
    s_beta.set_boundary(tk::spline::first_deriv,0.0,tk::spline::first_deriv,0,false);
    s_alfa.set_boundary(tk::spline::first_deriv,0.0,tk::spline::first_deriv,0,false);

    // T needs to be sorted, strictly increasing
    s_x.set_points(t,cartpos_x);
    s_y.set_points(t,cartpos_y);
    s_z.set_points(t,cartpos_z);

    // T needs to be sorted, strictly increasing
    s_gammar.set_points(t,cartpos_gammar);
    s_beta.set_points(t,cartpos_beta);
    s_alfa.set_points(t,cartpos_alfa);


}

KDL::JntArray trajectory_generate::homing(KDL::JntArray jointpositions_s,
                     double t,double t_s,double t_f){


    KDL::JntArray jointpositions_standard(6);
    jointpositions_standard.data <<PI/2,-PI/2,PI/2,PI/2,PI/2,PI/2;
    KDL::JntArray jointpositions(6);

        for(unsigned int i=0;i<6;i++){

            jointpositions(i) = jointpositions_s(i)
             +(3/((t_f-t_s)*(t_f-t_s)))*(jointpositions_standard(i)-jointpositions_s(i))*(t-t_s)*(t-t_s)
             -(2/((t_f-t_s)*(t_f-t_s)*(t_f-t_s)))*(jointpositions_standard(i)-jointpositions_s(i))*(t-t_s)*(t-t_s)*(t-t_s);
        }

    return jointpositions;
}
