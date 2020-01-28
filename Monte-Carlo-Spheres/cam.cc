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

Cam::Cam(Vector3d eye, Vector3d look, Vector3d up, vector<double> bounds,double height,double width, double focalLength) :
eye(eye), look(look), up(up), bounds(bounds), height(height), width(width), focalLength(focalLength)
{
}

Ray Cam :: pixelRay(int i, int j){
    Wv = eye-look; Wv = Wv.normalized();
    Uv = up.cross(Wv);Uv = Uv.normalized() ;
    Vv = Wv.cross(Uv);
    double px = i/(width-1.0)* (bounds[1]-bounds[0])+ bounds[0];
    double py = j/(height-1.0) * (bounds[2]-bounds[3])+ bounds[3];
	//cout << "px " << px << " py " << py << "\n";
    Vector3d point;
    Vector3d direction;
    point = eye + (focalLength * Wv) + (px * Uv) + (py * Vv);
    direction = point - eye; direction = direction.normalized();
    Ray ray = Ray(point,direction);
	return ray;
	
}


    
    
    
    
    
    
    
    
    
