#ifndef OBJECT_H_INCLUDED
#define OBJECT_H_INCLUDED
#include<ostream>
#include <iostream>
#include<vector>
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Eigen"
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "./Eigen/eigen-eigen-323c052e1731/Eigen/Geometry"
#include <memory>
//#include "HitRecord.h"


//using namespace std;



struct Material {
    double alpha;
    double illum;
    Eigen::Vector3d albedo;
    Eigen::Vector3d emitance;
    Ray scatter;
    std::string name;
    int matType;
    
    Material& operator=(Material mat){
        albedo = mat.albedo;
        name = mat.name;
        matType = mat.matType;
        return *this;
    }
    
};

class Object;

struct HitRecord {
    bool weCollidedBoys;
    double t;
    Eigen::Vector3d hitpoint;
    std::shared_ptr<Object> bestObj;
    Eigen::Vector3d normal;
    
    HitRecord& operator=(HitRecord h1){
        weCollidedBoys=h1.weCollidedBoys;
        t = h1.t;
        hitpoint = h1.hitpoint;
        bestObj = h1.bestObj;
		normal= h1.normal;
        return *this;
    }
    
};

class Object {
public:
    Material material;
    virtual Eigen::Vector3d getNormal(Eigen::Vector3d)=0;
    virtual HitRecord getCollision(Ray ray)=0;
	
};




#endif /*OBJECT_H_INCLUDED*/
