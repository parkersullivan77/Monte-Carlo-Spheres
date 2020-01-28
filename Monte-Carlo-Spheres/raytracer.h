#ifndef RAYTRACER_H_INCLUDED
#define RAYTRACER_H_INCLUDED
#include<ostream>
#include <iostream>
#include "cam.h"
#include "sphere.h"
#include "light.h"
#include "triangle.h"
#include "object.h"
#include <memory>
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Eigen"


class RayTracer{
	public:
    
	//methods
	void read(std::string);
	Eigen::Vector3d raySphereRGB(Ray,Eigen::Vector3d,double);
    void render(std::string);
	HitRecord checkHit(Ray);
   	void makeTriangles(std::vector<std::string>);
   	Eigen::MatrixXd makeTranslate(std::vector<std::string>);
	Eigen::MatrixXd makeRotate(std::vector<std::string>);
   	Eigen::MatrixXd makeScale(std::vector<std::string>);
   	Eigen::MatrixXd makeVertices(std::string);
	void triangulo(Eigen::MatrixXd,std::string);
	Material getMaterial(std::vector<Material>,std::string);
    Eigen::Vector3d random_unit();
    
    //ClassVariables
	Cam cam;
	std::vector<Light> lights;
    double level;
	//std::vector<Sphere> spheres;
    Eigen::Vector3d ambient;
    std::vector<std::shared_ptr<Triangle>> triangles;
	std::vector<std::shared_ptr<Object>> objList;
    double cutOff;
    int samples;
    
	//Constructors
	RayTracer(std::string &);
    	RayTracer() = default;
};
#endif /*RAYTRACER_H_INCLUDED*/

