#include <iostream>
//
//#include <ccd.hpp>
#include <vector>
//#include<Utils.hpp>
//#include<subfunctions.h>

#include <array>
#include <tight_inclusion/inclusion_ccd.hpp>
#include <tight_inclusion/Rational.hpp>
#include <fstream>

using namespace inclusion_ccd;
void case_check(){
  

    const Eigen::Vector3d a0s(0.1, 0.1, 0.1);
    const Eigen::Vector3d a1s(0, 0, 1);
    const Eigen::Vector3d a0e(1, 0, 1);
    const Eigen::Vector3d a1e(0, 1, 1);
    const Eigen::Vector3d b0s(0.1, 0.1, 0.1);
    const Eigen::Vector3d b1s(0, 0, 0);
    const Eigen::Vector3d b0e(0, 1, 0);
    const Eigen::Vector3d b1e(1, 0, 0);
    

    bool res;
    std::array<double, 3> err={{-1,-1,-1}};
    double ms=1e-8;
    double toi;
    const double tolerance=1e-6;
    const double t_max=1;
    const int max_itr=1e6;
    double output_tolerance;
    const int CCD_TYPE=1;
    res=inclusion_ccd::edgeEdgeCCD_double(
        a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,
        err,ms,toi,tolerance,t_max,max_itr,output_tolerance,CCD_TYPE
    );
    if (!using_rational_method()){
        std::cout<<"the double ccd result is "<<res<<std::endl;
    }
    else{
    std::cout<<"the ccd using rational core result is "<<res<<std::endl;
    }
    }


int main(int argc, char* argv[])
{
//    inclusion_ccd::Rational a;
    case_check();
    

    std::cout<<"done"<<std::endl;

    return 1;
}
