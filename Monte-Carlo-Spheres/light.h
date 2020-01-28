#ifndef LIGHT_H_INCLUDED
#define LIGHT_H_INCLUDED
#include<ostream>
#include <iostream>
#include<vector>
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Eigen"

class Light{
	public:
    Eigen::Vector3d pos;
    Eigen::Vector3d Cv;    
    Light(Eigen::Vector3d,Eigen::Vector3d);
    Light() = default;
    ~Light() = default;
};
#endif /*LIGHT_H_INCLUDED*/
