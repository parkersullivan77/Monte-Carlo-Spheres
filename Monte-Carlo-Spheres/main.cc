#include "raytracer.h"
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

using namespace std;

int main(int argc, char* argv[]){
	 try{
          string file = argv[1];
          cout << "hello\n";
          RayTracer r = RayTracer(file);
          cout << "chicken\n";
          r.render(argv[2]);
    }
     catch (string err) {
            cerr << "ERROR: " << err << '\n';
        }
	return 0;
};
