#include "triangle.h"
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

Triangle::Triangle(Vector3d xV,Vector3d yV,Vector3d zV) : xV(xV),yV(yV),zV(zV)
{
		//norm = getNormal(xV);
		//tresNorms[0] = norm;
		//tresNorms[1]= norm;
		//tresNorms[2]= norm;
		
}
HitRecord Triangle :: getCollision(Ray ray){		
	Vector3d firstCol;
	Vector3d secCol;
	Vector3d thrdCol;
	Vector3d swap;
	double detM;
	double detM1;
	double beta;
	HitRecord h1;
	h1.weCollidedBoys = false;
	ray.direction = ray.direction.normalized();
	swap = xV-ray.point;
	firstCol =xV-yV;
	secCol = xV-zV;
	thrdCol = ray.direction;
	MatrixXd MM(3,3);
	MatrixXd MMs1(3,3);
	MM << firstCol,secCol,thrdCol;
	MMs1=MM;
	MMs1.col(0) = swap;
	// compute beta (and necessary determinantes)
	detM = MM.determinant();
	detM1= MMs1.determinant();
	beta = detM1/detM;
	//cout << "BETA" << beta << endl;
	// check if beta is not ok and exit if so
	if(beta >= 0){
		// compute gamma (and necessary determinant)
		double gamma;
		double detM2;
		MatrixXd MMs2(3,3);
		MMs2=MM;
		MMs2.col(1) = swap;
		detM2 = MMs2.determinant();
		gamma = detM2/detM;
		//cout << "GAMMA" << gamma << endl;
		if(gamma >=0 && beta+gamma <=1){
		// check gamma and combination and exit if not ok
		// compute t and necessary determinant
			double stValue;
			double detM3;
			MatrixXd MMs3(3,3);
			MMs3=MM;
			MMs3.col(2) = swap;
			detM3 = MMs3.determinant();
			stValue = detM3/detM;
			if(stValue > 0.000001){
				h1.normal = ((1-beta-gamma)*tresNorms[0] + beta*tresNorms[1] + gamma*tresNorms[2]).normalized();
                if(h1.normal.dot(ray.direction) > 0) {
                    h1.normal = -h1.normal;   
                }
                h1.t = stValue;
                h1.hitpoint = ray.point + stValue * ray.direction;
                h1.weCollidedBoys = true;
			}
			//cout << "stValue" << stValue << endl;
		}
	}
    return h1;
        
}

Vector3d Triangle :: getNormal(Vector3d point){
		Vector3d normal;
		normal = (xV-yV).cross(xV-zV);
		normal = normal.normalized();
		return normal;
}
