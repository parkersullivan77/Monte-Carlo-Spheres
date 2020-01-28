#ifndef CAM_H_INCLUDED
#define CAM_H_INCLUDED
#include "ray.h"
#include "sphere.h"
#include<ostream>
#include <iostream>
#include<vector>
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Eigen"


class Cam{
	public:
 
	//Class Variables
	Eigen::Vector3d eye;
    Eigen::Vector3d look;
    Eigen::Vector3d up;
    Eigen::Vector3d Wv;
    Eigen::Vector3d Uv;
    Eigen::Vector3d Vv;
    Eigen::Vector3d Ev;
    std::vector<double> bounds;
    double height;
    double width;
    double focalLength;

     
    //methods
    Ray pixelRay(int,int);
    
    //Constructors
    Cam(Eigen::Vector3d,Eigen::Vector3d,Eigen::Vector3d,std::vector<double>,double,double,double);
	Cam() = default;
    ~Cam() = default;
	
 
   
};
#endif /*CAM_H_INCLUDED*/
