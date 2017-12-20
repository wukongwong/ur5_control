#include <Eigen/Dense>
#include <kdl/frames_io.hpp>
#include <kdl/kdl.hpp>
#include <kdl/jntarray.hpp>

typedef Eigen::Matrix<double,3,3> Matrix33d;
typedef Eigen::Matrix<double,3,1> Matrix31d;
typedef Eigen::Matrix<double,4,1> Matrix41d;

// computes the distance between two cubes iCGAL::to_doublen R^3 using double
// as input type and some internal EXACT floating point type
#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Polytope_distance_d_traits_3.h>
#include <CGAL/Homogeneous.h>

#include <iostream>
#include <cassert>
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif
// use an EXACT kernel...
typedef CGAL::Homogeneous<ET>                 K;
typedef K::Point_3                            Point;
// ...and the traits class based on the exact kernel
typedef CGAL::Polytope_distance_d_traits_3<K> Traitss;
typedef CGAL::Polytope_distance_d<Traitss>     Polytope_distance;




//pair1
/**
 * \brief this function is caculate the point when the two rigid body's diatance is nearest
 * @param cartpos - the end cartpos of ur
 * @param d - the diatance between the two rigid body
 * @return the point coordinate of the nearest point
 */
std::vector<Matrix31d> od(KDL::Frame cartpos,double &d)
{

  Matrix33d R;

  R<< cartpos.M.data[0],cartpos.M.data[1],cartpos.M.data[2],
      cartpos.M.data[3],cartpos.M.data[4],cartpos.M.data[5],
      cartpos.M.data[6],cartpos.M.data[7],cartpos.M.data[8];

  Matrix31d op0;
  op0<<cartpos.p(0),cartpos.p(1),cartpos.p(2);

//p1
  Matrix31d p0p1;
  Matrix31d p1_m;
  p0p1<<0,0.040,0.040;
  p1_m= op0+R*p0p1;

//p2
  Matrix31d p0p2;
  Matrix31d p2_m;
  p0p2<<0,0.040,-0.040;
  p2_m= op0+R*p0p2;

  //p3
  Matrix31d p0p3;
  Matrix31d p3_m;
  p0p3<<0,-0.040,0.040;
  p3_m= op0+R*p0p3;

//p4
  Matrix31d p0p4;
  Matrix31d p4_m;
  p0p4<<0,-0.040,-0.040;
  p4_m= op0+R*p0p4;

  //p5
  Matrix31d p0p5;
  Matrix31d p5_m;
  p0p5<<-0.12,0.040,0.040;
  p5_m= op0+R*p0p5;

  //p6
  Matrix31d p0p6;
  Matrix31d p6_m;
  p0p6<<-0.12,0.040,-0.040;
  p6_m= op0+R*p0p6;

  //p7
  Matrix31d p0p7;
  Matrix31d p7_m;
  p0p7<<-0.120,-0.040,0.040;
  p7_m= op0+R*p0p7;

  //p8
  Matrix31d p0p8;
  Matrix31d p8_m;
  p0p8<<-0.12,-0.040,-0.040;
  p8_m= op0+R*p0p8;

  std::ofstream ofile;
  // output the point
  ofile.open("/home/w/catkin_ws/devel/lib/publisher/point.txt",std::ios_base::app);
  ofile<< p1_m(0)<<" "<<p1_m(1)<<" "<<p1_m(2)<<std::endl
       << p2_m(0)<<" "<<p2_m(1)<<" "<<p2_m(2)<<std::endl
       << p3_m(0)<<" "<<p3_m(1)<<" "<<p3_m(2)<<std::endl
       << p4_m(0)<<" "<<p4_m(1)<<" "<<p4_m(2)<<std::endl
       << p5_m(0)<<" "<<p5_m(1)<<" "<<p5_m(2)<<std::endl
       << p6_m(0)<<" "<<p6_m(1)<<" "<<p6_m(2)<<std::endl
       << p7_m(0)<<" "<<p7_m(1)<<" "<<p7_m(2)<<std::endl
       << p8_m(0)<<" "<<p8_m(1)<<" "<<p8_m(2)<<std::endl;

  ofile.close();

  // the cube [0,1]^3
  Point P[8] = {Point(p4_m(0),p4_m(1),p4_m(2)),Point(p8_m(0),p8_m(1),p8_m(2)),Point(p3_m(0),p3_m(1),p3_m(2)),Point(p7_m(0),p7_m(1),p7_m(2)),
                Point(p2_m(0),p2_m(1),p2_m(2)),Point(p6_m(0),p6_m(1),p6_m(2)),Point(p1_m(0),p1_m(1),p1_m(2)),Point(p5_m(0),p5_m(1),p5_m(2))};

//  Point Q[8] = { Point(0.2,0.3,0), Point(0.2,0.3,0.2),Point(0.2,0.5,0), Point(0.2,0.5,0.2),
//                 Point(0.4,0.3,0), Point(0.4,0.3,0.2),Point(0.4,0.5,0), Point(0.4,0.5,0.2)};

  Point Q[8] = { Point(0.3,0.6,0.2), Point(0.3,0.6,0.4),Point(0.3,0.5,0.2), Point(0.3,0.5,0.4),
                 Point(0.4,0.6,0.2), Point(0.4,0.6,0.4),Point(0.4,0.5,0.2), Point(0.4,0.5,0.4)};;

  Polytope_distance pd(P, P+8, Q, Q+8);
  assert (pd.is_valid());
  // get squared distance (2,2,2)-(1,1,1))^2 = 3
  std::cout << "Squared distance: " <<
    CGAL::to_double (pd.squared_distance_numerator()) /
    CGAL::to_double (pd.squared_distance_denominator()) << std::endl;

  d =  CGAL::to_double (pd.squared_distance_numerator()) /
      CGAL::to_double (pd.squared_distance_denominator());
  d = sqrt(d);

 //output the desired cartendpos

  ofile.open("distance.txt",std::ios_base::app);
  ofile<< d
      <<std::endl;
  ofile.close();

  // get points that realize the distance
  Polytope_distance::Coordinate_iterator  coord_it;
//  std::cout << "p:"; // homogeneous point from first cube, (1,1,1,1)

  Matrix31d p;
  Matrix41d pp;

//  std::vector<double> p(3,0);

  int i=0;
  for (coord_it = pd.realizing_point_p_coordinates_begin();
       coord_it != pd.realizing_point_p_coordinates_end();
       ++coord_it){

//    std::cout << " " << *coord_it;
    pp(i) =CGAL::to_double(*coord_it);
    i++;

  }

  p<<pp(0)*(1/pp(3)),pp(1)*(1/pp(3)),pp(2)*(1/pp(3));

  ofile.open("p.txt",std::ios_base::app);
  ofile<< p(0)
      <<std::endl;
  ofile.close();

  std::cout << std::endl;

//  std::cout << "q:"; // homogeneous point from second cube, (2,2,2,1)
  Matrix31d q,pq;
  Matrix41d qq;
  int j=0;
  for (coord_it = pd.realizing_point_q_coordinates_begin();
       coord_it != pd.realizing_point_q_coordinates_end();
       ++coord_it){

//    std::cout << " " << *coord_it;
    qq(j) =CGAL::to_double(*coord_it);
    j++;
  }
  q<<qq(0)*(1/qq(3)),qq(1)*(1/qq(3)),qq(2)*(1/qq(3));

  std::cout << std::endl;
  std::vector<Matrix31d> geometric_info(2);
  geometric_info[0] = p;
  geometric_info[1] = q;
  double ddddd;
  pq = p-q;
  ddddd =sqrt( pq(0)*pq(0)+pq(1)*pq(1)+pq(2)*pq(2));
  std::cout << std::endl;
  return geometric_info;
}

//pair2

std::vector<Matrix31d> od1(KDL::Frame cartpos,double &d)
{

  Matrix33d R;

  R<< cartpos.M.data[0],cartpos.M.data[1],cartpos.M.data[2],
      cartpos.M.data[3],cartpos.M.data[4],cartpos.M.data[5],
      cartpos.M.data[6],cartpos.M.data[7],cartpos.M.data[8];

  Matrix31d op0;
  op0<<cartpos.p(0),cartpos.p(1),cartpos.p(2);

//p1
  Matrix31d p0p1;
  Matrix31d p1_m;
  p0p1<<0.02,0.04,-0.04;
  p1_m= op0+R*p0p1;

//p2
  Matrix31d p0p2;
  Matrix31d p2_m;
  p0p2<<0.02,0.04,0.04;
  p2_m= op0+R*p0p2;

  //p3
  Matrix31d p0p3;
  Matrix31d p3_m;
  p0p3<<0.02,-0.04,-0.04;
  p3_m= op0+R*p0p3;

//p4
  Matrix31d p0p4;
  Matrix31d p4_m;
  p0p4<<0.02,-0.04,0.04;
  p4_m= op0+R*p0p4;

  //p5
  Matrix31d p0p5;
  Matrix31d p5_m;
  p0p5<<-0.1,0.04,-0.04;
  p5_m= op0+R*p0p5;

  //p6
  Matrix31d p0p6;
  Matrix31d p6_m;
  p0p6<<-0.1,0.04,0.04;
  p6_m= op0+R*p0p6;

  //p7
  Matrix31d p0p7;
  Matrix31d p7_m;
  p0p7<<-0.1,-0.04,-0.04;
  p7_m= op0+R*p0p7;

  //p8
  Matrix31d p0p8;
  Matrix31d p8_m;
  p0p8<<-0.1,-0.04,0.04;
  p8_m= op0+R*p0p8;

  Matrix31d p86;

  p86 = p8_m-p6_m;
  double l= p86.norm();
  std::cout<<"L"<<l<<std::endl;

  std::ofstream ofile;

  // output the point
  ofile.open("/home/w/catkin_ws/devel/lib/publisher/point.txt",std::ios_base::app);
  ofile<< p1_m(0)<<" "<<p1_m(1)<<" "<<p1_m(2)<<std::endl
          << p2_m(0)<<" "<<p2_m(1)<<" "<<p2_m(2)<<std::endl
             << p3_m(0)<<" "<<p3_m(1)<<" "<<p3_m(2)<<std::endl
                << p4_m(0)<<" "<<p4_m(1)<<" "<<p4_m(2)<<std::endl
                   << p5_m(0)<<" "<<p5_m(1)<<" "<<p5_m(2)<<std::endl
                      << p6_m(0)<<" "<<p6_m(1)<<" "<<p6_m(2)<<std::endl
                         << p7_m(0)<<" "<<p7_m(1)<<" "<<p7_m(2)<<std::endl
                            << p8_m(0)<<" "<<p8_m(1)<<" "<<p8_m(2)<<std::endl;


  ofile.close();


  // the cube [0,1]^3

  Point P[8] = {Point(p4_m(0),p4_m(1),p4_m(2)),Point(p8_m(0),p8_m(1),p8_m(2)),Point(p3_m(0),p3_m(1),p3_m(2)),Point(p7_m(0),p7_m(1),p7_m(2)),
                Point(p2_m(0),p2_m(1),p2_m(2)),Point(p6_m(0),p6_m(1),p6_m(2)),Point(p1_m(0),p1_m(1),p1_m(2)),Point(p5_m(0),p5_m(1),p5_m(2))};



//  Point Q[8] = { Point(0.25,0.65,0.1), Point(0.25,0.65,0.3),Point(0.25,0.55,0.1), Point(0.25,0.55,0.3),
//                 Point(0.35,0.65,0.1), Point(0.35,0.65,0.3),Point(0.35,0.55,0.1), Point(0.35,0.55,0.3)};

  Point Q[8] = { Point(0.3,0.6,0.2), Point(0.3,0.6,0.4),Point(0.3,0.5,0.2), Point(0.3,0.5,0.4),
                 Point(0.4,0.6,0.2), Point(0.4,0.6,0.4),Point(0.4,0.5,0.2), Point(0.4,0.5,0.4)};;

  Polytope_distance pd(P, P+8, Q, Q+8);
  assert (pd.is_valid());
  // get squared distance (2,2,2)-(1,1,1))^2 = 3
  std::cout << "Squared distance: " <<
    CGAL::to_double (pd.squared_distance_numerator()) /
    CGAL::to_double (pd.squared_distance_denominator()) << std::endl;

  d =  CGAL::to_double (pd.squared_distance_numerator()) /
      CGAL::to_double (pd.squared_distance_denominator());
  d = sqrt(d);

//  Point pp,qq;
//  pp =pd.realizing_point_p();
//  qq = pd.realizing_point_q();
//  std::cout<<"pp"<<pp<<std::endl;


 //output the desired cartendpos

  ofile.open("distance.txt",std::ios_base::app);
  ofile<< d
      <<std::endl;
  ofile.close();

  // get points that realize the distance
  Polytope_distance::Coordinate_iterator  coord_it;
  std::cout << "p:"; // homogeneous point from first cube, (1,1,1,1)

  Matrix31d p;
  Matrix41d pp;


  int i=0;
  for (coord_it = pd.realizing_point_p_coordinates_begin();
       coord_it != pd.realizing_point_p_coordinates_end();
       ++coord_it){

    std::cout << " " << *coord_it;
    pp(i) =CGAL::to_double(*coord_it);
    i++;

  }

  p<<pp(0)*(1/pp(3)),pp(1)*(1/pp(3)),pp(2)*(1/pp(3));

  ofile.open("p.txt",std::ios_base::app);
  ofile<< p(0)
      <<std::endl;
  ofile.close();

  std::cout << std::endl;

  std::cout << "q:"; // homogeneous point from second cube, (2,2,2,1)
  Matrix31d q,pq;
  Matrix41d qq;
  int j=0;
  for (coord_it = pd.realizing_point_q_coordinates_begin();
       coord_it != pd.realizing_point_q_coordinates_end();
       ++coord_it){

    std::cout << " " << *coord_it;
    qq(j) =CGAL::to_double(*coord_it);
    j++;

  }
  q<<qq(0)*(1/qq(3)),qq(1)*(1/qq(3)),qq(2)*(1/qq(3));


  std::cout << std::endl;

  std::vector<Matrix31d> geometric_info(2);
  geometric_info[0] = p;
  geometric_info[1] = q;

  return geometric_info;
}



