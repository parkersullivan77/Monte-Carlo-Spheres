#ifndef RAY_H_INCLUDED
#define RAY_H_INCLUDED
#include<ostream>
#include <iostream>
#include<vector>
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Eigen"

class Ray{
	public:
    Eigen::Vector3d point;
    Eigen::Vector3d direction;    
    Ray(Eigen::Vector3d,Eigen::Vector3d);
    Ray() = default;
    ~Ray() = default;
};
#endif /*RAY_H_INCLUDED*/
