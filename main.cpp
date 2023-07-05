#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <chrono>
#include "Vector.cpp"

#include "lbfgs.c"

#include <sstream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define VOLUME_AIR 0.7
#define VOLUME_FLUID 0.3

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
    double area() const {
        //slide 9 
        if (vertices.size()<3) return 0; // not a polygon (not enough vertices)
        double result = 0;
        for (int i=0; i< vertices.size(); i++){
            //Vector A = vertices[i];
            //Vector B = vertices[(i+1)%vertices.size()];
            result += vertices[i][0]*vertices[(i+1)%vertices.size()][1] - vertices[i][1]*vertices[(i+1)%vertices.size()][0];
        }
        return std::abs(result/2);
    }

    double integrate_squared_distance(const Vector& P) const {
        if (vertices.size()<3) return 0; // not a polygon (not enough vertices)
        //slide 34
        double value =0; 
        for(int i = 1; i<vertices.size()-1;i++){
            Vector triangle[3] = {vertices[0], vertices[i], vertices[i+1]};
            
            double local_value = 0; 
            for (int k = 0; k<3; k++){
                for (int l =k; l<3; l++){
                    local_value += dot(triangle[k]-P, triangle[l]-P);
                }
            }

            Vector e1= triangle[1]-triangle[0] ; 
            Vector e2= triangle[2]-triangle[0];

            double area_triangle = 0.5 * abs(e1[1]*e2[0] - e1[0]*e2[1]);
            value += local_value*area_triangle/6.;
        }
        return value; 
    }

    Vector centroid() const {
        Vector c(0,0,0);
        int N = vertices.size();
        for (int i=0; i< vertices.size(); i++){
            c[0] += (vertices[i][0] + vertices[(i + 1) % N][0])* (vertices[i][0] * vertices[(i + 1) % N][1] - vertices[(i + 1) % N][0] * vertices[i][1]);
            c[1] += (vertices[i][1] + vertices[(i + 1) % N][1])* (vertices[i][0] * vertices[(i + 1) % N][1] - vertices[(i + 1) % N][0] * vertices[i][1]);
        }
        return c / (6 * area());
    }

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


int sgn(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0) {
    int W = 500, H = 500;
    std::vector<unsigned char> image(W*H * 3, 255);
// #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W-1., std::max(0., W * bminx));
        bminy = std::min(H-1., std::max(0., H * bminy));
        bmaxx = std::max(W-1., std::max(0., W * bmaxx));
        bmaxy = std::max(H-1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
                    int sign = sgn(det);
                    if (prevSign == 0) prevSign = sign; else
                        if (sign == 0) sign = prevSign; else
                        if (sign != prevSign) {
                            isInside = false;
                            break;
                        }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
                    double distEdge = std::abs(det)/ edgeLen;
                    double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
                    if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) { // want the cells blue
                //   // the N first particles may represent fluid, displayed in blue
                //      image[((H - y - 1)*W + x) * 3] = 0;
                //      image[((H - y - 1)*W + x) * 3 + 1] = 0;
                //      image[((H - y - 1)*W + x) * 3 + 2] = 255;
                
                    if (mindistEdge <= 2) {
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        image[((H - y - 1)*W + x) * 3 + 2] = 0;
                    }
                }
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

class PowerDiagram{
public: 
    PowerDiagram(){
        const int N_disk = 50;
        disk.vertices.resize(N_disk);

        for (int i=0; i<N_disk; i++){
            double t = i/(double)N_disk * M_PI * 2;
            disk.vertices[i][0] = cos(t);
            disk.vertices[i][1] = sin(t);
            disk.vertices[i][2] = 0;
        }
    };

    Polygon disk; 
    PowerDiagram(const std::vector<Vector>& pts, const std::vector<double> &weights){
        points = pts; 
        this->weights = weights;
        const int N_disk = 50; 
        disk.vertices.resize(N_disk);

        for (int i=0; i< N_disk; i++){
            double t = i / (double)N_disk*M_PI;
            disk.vertices[i][0] = cos(t);
            disk.vertices[i][1] = sin(t);
            disk.vertices[i][2] = 0;
        }
    }

    Polygon clip_polygon_by_bissector(const Polygon& poly,int index_0, int index_i, const Vector& P0, const Vector& Pi) const{
        //Sutherland Hodgman algorithm 
        // to check if its inside outside or in between 

        Vector M = (P0+Pi)*0.5;
        Vector Mprime = M + ((weights[index_0]-weights[index_i])/(2.*(P0-Pi).norm2()))*(Pi-P0);
        Polygon result;
        result.vertices.reserve(poly.vertices.size() + 1);
        for (int i=0; i< poly.vertices.size(); i++){

            const Vector &A = (i==0)? poly.vertices[poly.vertices.size()-1]:poly.vertices[i-1]; 
            const Vector &B = poly.vertices[i];
     
            double t = dot(Mprime-A, Pi-P0)/dot(B-A, Pi-P0); 
            Vector P = A+ t*(B-A); // Point of intersection?

            if ((B - P0).norm2() -weights[index_0] < (B - Pi).norm2()-weights[index_i]){ // B is inside 
                if ((A - P0).norm2() -weights[index_0] > (A - Pi).norm2()-weights[index_i]){ // A is outside 
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else if ((A - P0).norm2() -weights[index_0] < (A - Pi).norm2()-weights[index_i]){ // A is inside 
               // if ((B - P0).norm2()-weights[index_0]  > (B - Pi).norm2()-weights[index_i]){ // B is outside 
                result.vertices.push_back(P);
                //}
            }
        }
        return result;
    }

    Polygon clip_polygon_by_edge(const Polygon &poly, const Vector& u, const Vector& v) const {
        
        Polygon result;
        result.vertices.reserve(poly.vertices.size() + 1);
        Vector N(v[1]-u[1], -v[0]+u[0], 0);

        for (int i=0; i<poly.vertices.size(); i++) {
            const Vector &A = (i==0)? poly.vertices[poly.vertices.size() - 1]:poly.vertices[i-1];
            const Vector &B = poly.vertices[i];
            double t = dot(u-A, N) / dot(B-A, N);
            Vector P = A + t*(B-A);

            if (dot(B-u, N) < 0) { //B is inside
                if (dot(A-u, N) > 0) { //A outside
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else if (dot(A-u, N) < 0) { //A is inside
                result.vertices.push_back(P);
            }
        }
        return result;
    }

    Polygon intersect_with_disk(const Polygon& polygon, const Vector& center, double radius) const {
        Polygon result(polygon);
        for (int i=0; i<disk.vertices.size(); i++) {
            const Vector& u = disk.vertices[i]*radius + center; 
            const Vector& v = disk.vertices[(i+1)%disk.vertices.size()]* radius + center;
            result = clip_polygon_by_edge(result, u, v);
        }
        return result;
    }

    Polygon compute_powerdiagram_cell(int idx){
        
        Polygon result; 
        result.vertices.resize(4);
        //powerdiagram restricted to the unit square between 0 and 1 
        //clockwise direction
        result.vertices[0] = Vector(0,0,0);
        result.vertices[1] = Vector(1,0,0);
        result.vertices[2] = Vector(1,1,0);
        result.vertices[3] = Vector(0,1,0);

        // we want to clip it by bisector
        for (int i=0; i< points.size(); i++){
            if (i==idx) continue;
            result = clip_polygon_by_bissector(result,idx, i, points[idx], points[i]);
        }

        result = intersect_with_disk(result,points[idx], sqrt(weights[idx] - weights[weights.size()-1]));

        return result;
    }
    void compute(){ // compute polygon 
        powerdiagram.resize(points.size());
        for (int i=0; i< points.size(); i++){
            powerdiagram[i] = compute_powerdiagram_cell(i);
        }    
    }
    void save(std::string filename){
        save_svg(powerdiagram,filename, "blue" );
    }
    std::vector<Polygon> powerdiagram;
    std::vector<Vector> points;
    std::vector<double> weights;

};

class OT{
public:

    //OT(std::vector<Vector> &pts, const std::vector<double> &lambdas): pts(pts), lambdas(lambdas) {};
    OT(){};
    OT(const std::vector<Vector> &pts, const std::vector<double> &lambdas){
        this->pts = pts;
        this->lambdas=lambdas;
    };

    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<OT*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        lbfgsfloatval_t fx = 0.0;
     
        for (int i = 0;i < n-1;i++) {
            solution.weights[i] = x[i];
        }
        solution.compute();

        double s1 =0;
        double s2 =0;
        double s3 =0;
        double estimated_volume_fluid = 0;

        for (int i = 0;i < n+1;i++) {
            // slide 32 its the function g (W)
            //third term
            double cell_area = solution.powerdiagram[i].area();
            g[i]= -(lambdas[i] - cell_area);
            s3 = s3 + lambdas[i]*x[i];
            //second term 
            s2 = s2 - (x[i]*cell_area);
            //first term
            s1 = s1 + solution.powerdiagram[i].integrate_squared_distance(solution.points[i]);
            estimated_volume_fluid += cell_area;

        }
        
        fx  = s1+ s2 + s3;

        double estimated_volume_air = 1. - estimated_volume_fluid;

        //Rest of formula for fluid 

        fx += x[n+1]*(VOLUME_AIR - estimated_volume_air);
        g[n+1]= -(VOLUME_AIR - estimated_volume_air);

        return -fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<OT*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        for (int i = 0;i < n-1;i++) {
            solution.weights[i] = x[i];
        }
        solution.compute();

        double max_diff=0;
        for(int i=0; i< n-1; i++){
            double current_area = solution.powerdiagram[i].area();
            double desired_area = lambdas[i];
            max_diff = std::max(max_diff, std::abs(current_area-desired_area));
        }

        std::cout<<"fx:"<< fx<<"\tmax difference = "<< max_diff<<"\t gnorm"<< gnorm <<std::endl;

        return 0;
    }


    void solve(){
        //optimal transport problem solution 
        solution.points = pts;
        solution.weights.resize(pts.size()+1); 
        std::fill(solution.weights.begin(), solution.weights.end(), 1.0);
        solution.weights[solution.weights.size()-1] = 0.7; 

        double fx =0;
        // LBFGS... copied from sample.cpp 
        int ret = lbfgs(pts.size()+1, &solution.weights[0], &fx, _evaluate, _progress, this, NULL);
        solution.compute();

    }

    std::vector<Vector> pts;
    std::vector<double> lambdas;
    PowerDiagram solution;


};

class Fluid {
public:
    Fluid(){};

    Fluid(int N){
        particles.resize(N);
        for (int i=0; i<N; i++) {
            particles[i] = Vector(rand()/(double)RAND_MAX, rand()/(double)RAND_MAX);
        }
        velocities.resize(N, Vector(0, 0, 0));
    }

    void stepFluid() {

        otsolver.pts = particles;
        otsolver.lambdas = std::vector<double>(particles.size(), 1./particles.size() * VOLUME_FLUID);
        otsolver.solve();

        const double mass_particle = 200;
        const double epsilon2 = 0.004*0.004;
        const double dt = 0.002;

        for (int i=0; i<particles.size(); i++) {
            Vector gravity = Vector(0, -9.81, 0) * mass_particle;
            Vector centroid = otsolver.solution.powerdiagram[i].centroid();
            Vector otForce = 1. / epsilon2 * (centroid - particles[i]);
            Vector forces = gravity + otForce;
            velocities[i] += dt / mass_particle * forces;
            particles[i] += dt * velocities[i];
        }
    }

    void runFluid() {
        for (int i=0; i < 1000; i++) {
            stepFluid();
            save_frame(otsolver.solution.powerdiagram, "animation", i);
        }
    }

    OT otsolver;
    std::vector<Vector> points;
    std::vector<Vector> velocities;
    std::vector<Vector> particles;
    std::vector<double> weights;
};


int main() {
    Fluid fluid(10);
    fluid.runFluid();
    exit(0);
    std::vector<Vector> points(32); 
    std::vector<double> lambdas(32); 
    //powerdiagram1 - 30 
    //powerdiagram2 - 1024
    //powerdiagram3 - 256
    //powerdiagram4 - debut td 7 

    //std::vector<Vector> points(32); 
    //std::vector<double> lambdas(32); 
    for (int i=0; i< points.size(); i++){
        points[i][0] = rand()/ (double)RAND_MAX;
        points[i][1] = rand()/ (double)RAND_MAX;
        points[i][2] = 0;
        //weights[i] = 1; //powerdiagram4
        //weights[i] = rand()/ (double)RAND_MAX; //powerdiagram5
        lambdas[i] = 1./ points.size(); //powerdiagram6 with size 32 
    }
    OT ot(points, lambdas);
    auto start = std::chrono::high_resolution_clock::now();
    ot.solve();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout<<"duration = "<< duration.count() <<std::endl;
    ot.solution.save("powerdiagram5.svg");
    return 0;
}

// TD7:
// 1. Power diagrams:
// -> adjusting formulas and adding weights P and inside 
// 2. libLBFGS library moodle page for optimization
// 3. Optimal transport probleme:
// -> add libLBFGS library to the project
// -> implement evaluate function that computer the objetive function + gradient
 

 //TD8