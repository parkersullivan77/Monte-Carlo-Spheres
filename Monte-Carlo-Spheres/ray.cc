#include "ray.h"
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

Ray::Ray(Vector3d point, Vector3d direction): point(point),direction(direction)
{    
}

