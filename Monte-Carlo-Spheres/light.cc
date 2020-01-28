#include "light.h"
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

Light::Light(Vector3d pos, Vector3d Cv): pos(pos),Cv(Cv)
{    
}

