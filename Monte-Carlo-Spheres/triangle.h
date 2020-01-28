#ifndef TRIANGLE_H_INCLUDED
#define TRIANGLE_H_INCLUDED
#include<ostream>
#include <iostream>
#include<vector>
#include "ray.h"
#include "object.h"
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Eigen"
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Geometry"



class Triangle: public Object{
	public:
	//Class Variables
    Eigen::Vector3d xV;
    Eigen::Vector3d yV;
    Eigen::Vector3d zV;
	Eigen::Vector3d norm;
	std::vector<Eigen::Vector3d> tresNorms;
    Eigen::Vector3i vertIdx;
    int x,y,z;
    
	
	
    //Constructors
	Triangle(Eigen::Vector3d,Eigen::Vector3d,Eigen::Vector3d);
    Triangle() = default;
    ~Triangle() = default;
    
    //Mehtods
    HitRecord getCollision(Ray);
    Eigen::Vector3d getNormal(Eigen::Vector3d);
};
#endif /*TRIANGLE_H_INCLUDED*/
