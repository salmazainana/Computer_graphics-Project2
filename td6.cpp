#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <chrono>
#include "Vector.cpp"

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
    std::vector<Vector> vertices;
};  
 
// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
    FILE* f = fopen(filename.c_str(), "w+"); 
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i=0; i<polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \""); 
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke-width =\"5\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

class VoronoiDiagram{
public: 
    VoronoiDiagram(const std::vector<Vector>& pts){
        points = pts; 
    };

    Polygon clip_polygon_by_bissector(const Polygon& poly, const Vector& P0, const Vector& Pi){
        //Sutherland Hodgman algorithm 
        // to check if its inside outside or in between 
        Polygon result;
        for (int i=0; i< poly.vertices.size(); i++){

            Vector A = (i==0)? poly.vertices[poly.vertices.size()-1]:poly.vertices[i-1]; 
            const Vector &B = poly.vertices[i];
            // not the formula from standard algorithm but one on slide 61
            Vector M = (P0+Pi)*0.5;
            double t = dot(M-A, Pi-P0)/dot(B-A, Pi-P0); 
            Vector P = A+ t*(B-A); // Point of intersection?

            if ((B - P0).norm2() < (B - Pi).norm2()){ // B is inside 
                if ((A - P0).norm2() > (A - Pi).norm2()){ // A is outside 
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else if ((A - P0).norm2() < (A - Pi).norm2()){ // A is inside 
                if ((B - P0).norm2() > (B - Pi).norm2()){ // B is outside 
                    result.vertices.push_back(P);
                }
            }
        }
        return result;
    }

    Polygon compute_voronoi_cell(int idx){
        
        Polygon result; 
        result.vertices.resize(4);
        //voronoi restricted to the unit square between 0 and 1 
        //clockwise direction
        result.vertices[0] = Vector(0,0,0);
        result.vertices[1] = Vector(0,1,0);
        result.vertices[2] = Vector(1,1,0);
        result.vertices[3] = Vector(1,0,0);

        // we want to clip it by bisector
        for (int i=0; i< points.size(); i++){
            if (i==idx) continue;
            result = clip_polygon_by_bissector(result, points[idx], points[i]);
        }

        return result;
    }
    void compute(){ // compute polygon 
        voronoi.resize(points.size());
        for (int i=0; i< points.size(); i++){
            voronoi[i] = compute_voronoi_cell(i);
        }    
    }
    void save(std::string filename){
        save_svg(voronoi,filename, "blue" );
    }
    std::vector<Vector> points; 
    std::vector<Polygon> voronoi; 

};

int main(){
    //1 - 30 
    //2- 1024
    //3 - 256
    std::vector<Vector> points(256); 
    for (int i=0; i< points.size(); i++){
        points[i][0] = rand()/ (double)RAND_MAX;
        points[i][1] = rand()/ (double)RAND_MAX;
        points[i][2] = 0;
    }
    VoronoiDiagram voronoi(points);
    voronoi.compute();
    voronoi.save("voronoi3.svg");
    return 0;
}