#ifndef SPHERE_H_INCLUDED
#define SPHERE_H_INCLUDED
#include<ostream>
#include <iostream>
#include<vector>
#include "ray.h"
#include "sphere.h"
#include "object.h"
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Eigen"
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Geometry"


class Sphere: public Object {
	public:
	//Class Variables
    Eigen::Vector3d xyz;
    //Eigen::Vector3d Ka;
    //Eigen::Vector3d Kd;
    //Eigen::Vector3d Ks;
    //Eigen::Vector3d Kr; 
    //Material mat;
    double r;
    Eigen::Vector3d best_pt;
    
	
	
    //Constructors
    Sphere(Eigen::Vector3d,double);
    Sphere() = default;
    ~Sphere() = default;
    
    //Mehtods
   HitRecord getCollision(Ray);
   Eigen::Vector3d getNormal(Eigen::Vector3d);

};
#endif /*SPHERE_H_INCLUDED*/
