#include "sphere.h"
#include "cam.h"
#include <string>
#include<cctype>
#include<vector>
#include <iostream>
#include <stdlib.h>
#include<getopt.h>
#include <sstream>
#include<fstream>
#include<ostream>
#include<map>
#include <bits/stdc++.h> 
#include <boost/algorithm/string.hpp>
#include <cmath>


using namespace Eigen;
using namespace std;


Sphere::Sphere(Vector3d xyz,double r) : xyz(xyz),r(r)
{
}


HitRecord Sphere :: getCollision(Ray ray){
    //cout << "IN SPHERE" << endl;
	HitRecord h1;
        h1.weCollidedBoys = false;
		double v;
        double csq;
        double disc;
        double d;
        double t; 
        double pt;
        Vector3d Tv;
        Tv = xyz-ray.point;
        v = Tv.dot(ray.direction);
        csq =Tv.dot(Tv);
        disc = r*r - (csq -(v*v));
        if(disc < 0){
            return h1;
        }
        else
	    h1.weCollidedBoys = true;
            d = sqrt(disc);
            t = v-d;
	    h1.t=t;
	    h1.hitpoint = ray.point + t*ray.direction;
	    h1.normal = getNormal(h1.hitpoint);
        return h1;
        
    }
Vector3d Sphere :: getNormal(Vector3d point){
		return point - this->xyz;
}
