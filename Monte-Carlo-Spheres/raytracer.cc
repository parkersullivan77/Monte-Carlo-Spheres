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
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <memory>
#include <cmath>
#include <random>
#include "object.h"

using namespace std;
using namespace Eigen;


RayTracer :: RayTracer(string &filename){
    read(filename);
}

void RayTracer :: read(string filename){
    vector<double> bounds;
    double r;
    string driverLine;
    vector<string> splitLines;
    Vector3d eye;
    Vector3d look;
    Vector3d up;
    double height;
    double width;
    double focalLength;
    ifstream input(filename);
    while(getline(input,driverLine)){
    boost::split(splitLines,driverLine,boost::is_any_of(" "));   
    if( splitLines[0] == "eye"){
        eye << stod(splitLines[1]),stod(splitLines[2]),stod(splitLines[3]);
    
        
    }
    if(splitLines[0] == "look"){
        look << stod(splitLines[1]),stod(splitLines[2]),stod(splitLines[3]);
    }
    if(splitLines[0]== "up"){
        up << stod(splitLines[1]),stod(splitLines[2]),stod(splitLines[3]);
    }
    if(splitLines[0] == "bounds"){
        bounds.push_back(stod(splitLines[1]));
        bounds.push_back(stod(splitLines[2]));
        bounds.push_back(stod(splitLines[3]));
        bounds.push_back(stod(splitLines[4]));
    }
    if(splitLines[0] == "res"){
        height = stod(splitLines[1]);
        width = stod(splitLines[2]);
    }
    
    /*if(splitLines[0] == "light"){
        Vector3d pos(stod(splitLines[1]),stod(splitLines[2]),stod(splitLines[3]));
        Vector3d Cv(stod(splitLines[5]),stod(splitLines[6]),stod(splitLines[7]));
        lights.push_back(Light(pos,Cv));
    }*/
    
    if(splitLines[0] == "sphere"){
    
	Vector3d null(0,0,0);
    Vector3d xyz(stod(splitLines[1]),stod(splitLines[2]),stod(splitLines[3]));
    r = (stod(splitLines[4]));
    shared_ptr<Sphere> currentSphere (new Sphere(xyz,r));
    currentSphere->material.name = "Sphere";	
    if(stoi(splitLines[8]) == 1){
         currentSphere->material.matType = 1;
         Vector3d albedo(stod(splitLines[5]),stod(splitLines[6]),stod(splitLines[7]));
         currentSphere->material.albedo = albedo;
         currentSphere->material.emitance = null;
         //cout << "emittance\n" << currentSphere->material.emitance << "\n" << "albedo\n" <<currentSphere->material.albedo << "\n"; 
    }
    
    if(stoi(splitLines[8]) == 2){
        currentSphere->material.matType = 2;
        Vector3d albedo(stod(splitLines[5]),stod(splitLines[6]),stod(splitLines[7]));
        currentSphere->material.albedo = albedo;
        currentSphere->material.emitance = null;
        //cout << "emittance\n" << currentSphere->material.emitance << "\n" << "albedo\n" <<currentSphere->material.albedo << "\n"; 
    }
    
    if(stoi(splitLines[8]) == 0){
        Vector3d emitted(stod(splitLines[5]),stod(splitLines[6]),stod(splitLines[7]));
        currentSphere->material.matType = 0;
        currentSphere->material.emitance = emitted;
    }
        objList.push_back(currentSphere);
    }
	
    if(splitLines[0] == "ambient"){
        ambient << stod(splitLines[1]),stod(splitLines[2]),stod(splitLines[3]);
    }
    if(splitLines[0] == "d"){
        focalLength = stod(splitLines[1]);
    }
    //if(splitLines[0] == "model"){ 
    //    makeTriangles(splitLines);
    // }
    if(splitLines[0] == "recursionlevel"){
        level = stod(splitLines[1]);
    }
    if(splitLines[0] == "samples"){
        samples = stoi(splitLines[1]);
    }
}
    /*cout << "look\n" << look << "\n";
    cout << "eye\n" << eye << "\n";
    cout << "up\n" << up << "\n";
    cout << "Height Width " << height << " " << width << "\n"; 
    cout << "focalLength " << focalLength << "\n";
    cout << "recursionlevel " << level << "\n";*/
   cam = Cam(eye,look,up,bounds,height,width,focalLength);


}


/*void RayTracer :: makeTriangles(vector<string> modelSpecs){
    vector<string> tValues;
    MatrixXd tM; 
    MatrixXd rM;
    MatrixXd sM;
    MatrixXd MM;
    MatrixXd l;
    MatrixXd transformedVertices;
    MatrixXd unrotatedVertices;
    
    tValues.push_back(modelSpecs[1]);
    tValues.push_back(modelSpecs[2]);
    tValues.push_back(modelSpecs[3]);
    tValues.push_back(modelSpecs[4]);
    tValues.push_back(modelSpecs[5]);
    tValues.push_back(modelSpecs[6]);
    tValues.push_back(modelSpecs[7]);
    tValues.push_back(modelSpecs[8]);
    cutOff = stod(modelSpecs[9]);
    rM = makeRotate(tValues);
    tM = makeTranslate(tValues);
    sM = makeScale(tValues);
    MM = tM * sM * rM;
        
    string filename = modelSpecs[10];
    unrotatedVertices = makeVertices(filename);
    //cout << unrotatedVertices.col(0);
    //cout << MM.col(0);
    //cout << "MM\n" << MM << endl; 
    transformedVertices = MM * unrotatedVertices;
    //cout << "unrotatedVertices\n" << unrotatedVertices << endl;
    //cout << "vertices\n" << transformedVertices << endl;
    //triangles = triangulo(transformedVertices,filename);
    triangulo(transformedVertices,filename);
    
}

MatrixXd RayTracer::makeTranslate(vector<string> tValues){
        double tx = stod(tValues[5]);
        double ty = stod(tValues[6]);
        double tz = stod(tValues[7]);
        MatrixXd translate(4,4);
        translate << 1, 0, 0, tx,
                    0, 1, 0, ty,
                    0, 0, 1, tz,
                    0, 0, 0, 1;
        //cout << translate << '\n';
        return translate;
}

MatrixXd RayTracer::makeScale(vector<string> tValues){
        double scale = stod(tValues[4]);
        MatrixXd sM(4,4);
        sM << scale,0,0,0,
            0,scale,0,0,
            0,0,scale,0,
            0,0,0,1;
        //cout << sM;
        return sM;
}

MatrixXd RayTracer::makeRotate(vector<string> tValues){
        double wx = stod(tValues[0]);
        double wy = stod(tValues[1]);
        double wz = stod(tValues[2]);
        double theta = stod(tValues[3]);
        double arad = (theta/180)*M_PI;
        double ca = cos(arad);
        double sa = sin(arad);
        Vector3d Wv(wx,wy,wz);
        Wv = (1/Wv.norm())*Wv;
        Vector3d Mv;
        Mv = Wv;
        //cout << "Mv\n " << Mv << " \nWv\n" << Wv<<endl;
        ptrdiff_t i, j;
        double minOfM = Mv.minCoeff(&i,&j);
        //cout <<"i "<<i <<" j " << j << endl;
        //cout<< "MIN OF M " << Mv[i,j] << "\n";
        Mv[j,i] = 1;
        //cout << "EARLY MV\n" << Mv<<endl;
        Vector3d Uv;
        Uv = Wv.cross(Mv);
        Uv = (1/Uv.norm())*Uv;
        Vector3d Vv;
        Vv = Wv.cross(Uv);
        
        //cout << "MV" << '\n' << Mv << '\n' << "WV" << '\n' << Wv << '\n' << "Uv" << Uv << '\n' << "Vv" << Vv << '\n' ; 		
        
        MatrixXd rM(3,3);
        rM << Uv(0), Uv(1),Uv(2),
            Vv(0), Vv(1),Vv(2),
            Wv(0), Wv(1),Wv(2);
        //cout <<"RM"<< rM << '\n' << '\n';
        
        MatrixXd rZ(3,3);
        rZ << ca,-sa,0,
            sa,ca,0,
            0,0,1;
        //cout<< "RZ" <<rZ << '\n' << '\n';
        
        MatrixXd rMT(3,3);
        rMT = rM.transpose();
        //cout << "RMT" << rMT << '\n'<< '\n';
        
        MatrixXd R(3,3);
        R = rMT * rZ * rM;
        //cout << R <<'\n' << '\n';
        MatrixXd RHT(4,4);
            RHT << R(0,0),R(0,1),R(0,2),0,
                R(1,0),R(1,1),R(1,2),0,
                R(2,0),R(2,1),R(2,2),0,
                0,0,0,1;
            
        //cout << RHT << "\n";
    return RHT;
}


MatrixXd RayTracer::makeVertices(string modelFile){
    ifstream input(modelFile);
    string vertex;
    vector<string> vertices;
    vector<string> splitLines;
    vector<string> points;
    string mtlFile;
    if(input){
    while(getline(input,vertex)){
            boost::split(splitLines,vertex,boost::is_any_of(" "));
        if(splitLines[0] == "mtllib"){
                mtlFile = splitLines[1];
            }
            if(splitLines[0] == "v"){
                vertices.push_back(vertex);
            }
            
    }
}
    MatrixXd vM(4,vertices.size());
    //cout << vertices.size();
    for(int i =0; i < vertices.size(); i++){
    for(int j=0; j < 3; j++){
        boost::split(points,vertices[i],boost::is_any_of(" "));
        vM(j,i) = stod(points[j+1]);
    }
        vM(3,i) = 1;
    }
    //cout<< "VM\n" << vM.col(63) << "\n";
    return vM;
}

Material RayTracer :: getMaterial(vector<Material> materials,string matName){
    for(Material mat : materials){
        if(mat.name == matName){
            return mat;
        }
    }
}

void RayTracer::triangulo(MatrixXd transVert, string filename){
    vector<Triangle> tList;
    string mtlFile;
    string vertex;
    vector<string> splitLines;
    vector<string> splitFace;
    vector<string> getName;
    vector<string> t;
    vector<Material> materials;
    string mtlLines;
    vector<string> mtlSplit;
    // Read in the mtllib file, construct materials
    ifstream input(filename);
    if(input){
    while(getline(input,vertex)){
            boost::split(getName,vertex,boost::is_any_of(" "));
        if(getName[0] == "mtllib"){
                mtlFile = getName[1];
                }
        }
    }
    ifstream mtl (mtlFile);
    if(mtl){
        
        while(getline(mtl,mtlLines)){
            if(mtlLines.size() == 0) {
                continue;
            }
            boost::split(mtlSplit,mtlLines,boost::is_any_of(" "));
            
            if(mtlSplit[0] == "newmtl"){
                materials.push_back(*(new Material()));
                materials[materials.size()-1].name =  mtlSplit[1];
                //mat.name = mtlSplit[0];
            }
            if(mtlSplit[0] == "Ka"){
                materials[materials.size()-1].Ka << stod(mtlSplit[1]),stod(mtlSplit[2]),stod(mtlSplit[3]);
            }
            if(mtlSplit[0] == "Kd"){
                materials[materials.size()-1].Kd << stod(mtlSplit[1]),stod(mtlSplit[2]),stod(mtlSplit[3]);
            }
            if(mtlSplit[0] == "Ks"){
                materials[materials.size()-1].Ks << stod(mtlSplit[1]),stod(mtlSplit[2]),stod(mtlSplit[3]);
            }
            if(mtlSplit[0] == "illum"){
        if(stod(mtlSplit[1]) >= 3){
        materials[materials.size()-1].Kr =  materials[materials.size()-1].Ks;		
        }
        else{
        materials[materials.size()-1].Kr << 0,0,0;
        }
            
            }
            if(mtlSplit[0] == "Ns"){
                materials[materials.size()-1].alpha = stod(mtlSplit[1]);
            }
        }
    Material currentMaterial;
    vector<vector<shared_ptr<Triangle>>> vertexToTriangle;

    ifstream objLine(filename);
        if(objLine){
            //cout << "hello";
    while(getline(objLine,vertex)){
        boost::split(splitLines,vertex,boost::is_any_of(" "));
            //cout << "hello:" <<splitLines[0] << "\n";
            if(splitLines[0] == "usemtl"){
                currentMaterial = getMaterial(materials,splitLines[1]);// get the material associted with that name
                }
            if (splitLines[0]== "v"){
                vector<shared_ptr<Triangle>> t;
                vertexToTriangle.push_back(t);
            }
            if(splitLines[0] == "f"){
                string trimmedFace = vertex.substr(2);
                int x,y,z;
                vector<string> vert; 
                boost::split(splitFace,trimmedFace,boost::is_any_of(" "));
                for(string face : splitFace){
                    //cout << face << "\n";
                    boost::split(t,face,boost::is_any_of("/"));
                    vert.push_back(t[0]);
                }
                //cout << "hello" <<vert[0] << " " << vert[1] << " " << vert[2] << "\n";
                x = stoi(vert[0]);
                y = stoi(vert[1]);
                z = stoi(vert[2]);
                
                Vector3d xV;
                Vector3d yV;
                Vector3d zV;
                
                xV << transVert.col(x-1)(0),transVert.col(x-1)(1),transVert.col(x-1)(2);
                yV << transVert.col(y-1)(0),transVert.col(y-1)(1),transVert.col(y-1)(2);
                zV << transVert.col(z-1)(0),transVert.col(z-1)(1),transVert.col(z-1)(2);
                shared_ptr<Triangle> tri (new Triangle(xV,yV,zV));
                tri->vertIdx(0) =x-1;
                tri->vertIdx(1) =y-1 ;
                tri->vertIdx(2) =z-1;
                auto triangleNormal = tri->getNormal(Vector3d(0,0,0));
                tri->norm = triangleNormal;
                vertexToTriangle[x-1].push_back(tri);
                vertexToTriangle[y-1].push_back(tri);
                vertexToTriangle[z-1].push_back(tri);
                tri->material = currentMaterial;
                triangles.push_back(tri);
                objList.push_back(tri);
            
            }
            }
        }
        for(auto &face: triangles) {
            //face->tresNorms.push_back(face->norm);
            //face->tresNorms.push_back(face->norm);
            //face->tresNorms.push_back(face->norm);
            Vector3d sum = Vector3d(0,0,0);
            sum = face->norm;
            for(int i =0; i <3; i++){
                int index = face->vertIdx(i);
                for(auto otherTriangle : vertexToTriangle[index] ){
                        double angle = otherTriangle->norm.dot(face->norm);
                        if(angle > 1){
                            angle =1;
                        }
                        if(angle < -1){
                            angle = -1;
                        }
                        if(acos(angle)*180/M_PI < cutOff){
                            sum += otherTriangle->norm;
                        }
                        sum = sum.normalized();
                        face->tresNorms.push_back(sum);
                    }
            }
        
        }

        
    }
    
}*/
    
HitRecord RayTracer :: checkHit(Ray ray){
        //cout << "INHERE" << endl;
    shared_ptr<Object> closest;
        HitRecord h1;
        h1.t = 1000000000;
        h1.weCollidedBoys = false;

        
    for(shared_ptr<Object> obj : objList){
    HitRecord temp = obj->getCollision(ray);
    if(temp.weCollidedBoys && temp.t > 0){
            if(temp.t < h1.t){
                h1 = temp;
                h1.bestObj = obj;
                }
            } 
    }
        return h1;
}

Vector3d RayTracer :: random_unit(){
    while (true){
        double random1 = ((float)(rand())/ (float)(RAND_MAX))*2.0-1.0;
        double random2 =  ((float)(rand())/ (float)(RAND_MAX))*2.0-1.0;
        double random3 =  ((float)(rand())/ (float)(RAND_MAX))*2.0-1.0;
        //cout << random1 << " " <<random2 << " " << random3 << endl;
        Vector3d dir(random1,random2,random3);dir = dir.normalized();
        double p = dir.dot(dir);
        if(p < 1){
            return dir;
        }
    }
    
}

Vector3d RayTracer :: raySphereRGB(Ray ray,Vector3d currentAlbedo,double level){
    //cout << "hello" << endl;
    HitRecord h1 = checkHit(ray);
    shared_ptr<Object> closest;
    double hitp = h1.t;
    Vector3d null(0,0,0);
    double random4 = (float)(rand())/(float)(RAND_MAX);
        
        if(h1.weCollidedBoys){
            //cout << hitp << "\n";
            closest = h1.bestObj;
            Vector3d best_pt;
            best_pt = ray.point + hitp * ray.direction;
            Vector3d normal;
            int materialType = closest->material.matType;
            normal = closest->getNormal(best_pt); normal = normal.normalized();
            
            bool doesScatter = true;

            switch(materialType){
                case 0 : // light
                    doesScatter = false;
					break;
                case 1 : //lambertian
                {  
                    Vector3d dir;
                    dir = random_unit() + normal;
                    dir =dir.normalized();
					closest->material.scatter = Ray(best_pt,dir);
					break;
                }
                case 2 :  //Mirror
                {
                    doesScatter = true;
	                Vector3d Uinv;
					Vector3d refR;
                    
	                //Uinv = -1 * ray.direction;
					//refR = (2 * normal.dot(Uinv) * normal)- Uinv); refR = refR.normalized();
					refR = ray.direction - 2 * (ray.direction.dot(normal) * normal);
                    closest->material.scatter = Ray(best_pt,refR);
					break;
                    
                }
                
                default:
                    doesScatter = false;
                    return null.cwiseProduct(currentAlbedo);
            }
            //cout << "OG level" << level << "\n";
            //cout << currentAlbedo << "\n";
            //cout << doesScatter << "\n";
            if(level > 0 && doesScatter){
                currentAlbedo = currentAlbedo.cwiseProduct(closest->material.albedo);
                //cout << currentAlbedo << "\n";
                if(currentAlbedo.maxCoeff() < random4){
                    return null.cwiseProduct(currentAlbedo);
                }
                currentAlbedo = currentAlbedo/currentAlbedo.maxCoeff();
                level = level-1;
                Vector3d out = raySphereRGB(closest->material.scatter,currentAlbedo,level);
				return out;
            }
			else{
				return currentAlbedo.cwiseProduct(closest->material.emitance);
			}
        }
        else{
            return null.cwiseProduct(currentAlbedo);
        }
        
}


//parker is a big ole cutie -love, W
    


void RayTracer :: render(string filename){
    //cout << "INR\n";
    Vector3d image[(int)cam.width][(int)cam.height];
    Vector3d black;
    Vector3d currentAlbedo;
    Ray ray;
    black << 0,0,0;
    int numSamples =0;
    int maxSamples = samples;
    
    
    Vector3d res;
    ofstream out;
    out.open(filename);
    out << "P3\n";
    out << cam.width << " " << cam.height <<  " 255\n";
    //cout << "NEXT BREAK\n";
    for(int i = 0; i < cam.width; i++){
        for(int j = 0; j < cam.height;j++){
            image[i][j] = black;
        }
    }
    
    while(numSamples < (maxSamples+1)){ 
    
    cout << "DONE WITH " << numSamples << " SAMPLES\n";
    numSamples += 1;
    
      for(int i =0; i < cam.width; i++){
        for(int j = 0; j < cam.height; j++){
          ray = (cam.pixelRay(i,j));
          currentAlbedo<< 1.0,1.0,1.0;
          res = raySphereRGB(ray,currentAlbedo,level);
          if(res(0) != 0 || res(1) != 0 || res(2) != 0){
             for(int e = 0; e <= 2; e++){
               res(e) = sqrt(res(e));
               image[i][j](e) += res(e);
            } 
          }
        }
      }
    }
    
    for(int i =0; i< cam.width;i++){
        for(int j =0; j < cam.height;j++){
            for(int e =0; e<=2;e++){
            image[j][i](e) = max(0.0,min(255.0,round(255.0 * image[j][i](e)/numSamples)));
           }
            out << (int)image[j][i](0) << " " << (int)image[j][i](1) << " " << (int)image[j][i](2) << "\n"; 
        }
    }
    out.close();
}






    



