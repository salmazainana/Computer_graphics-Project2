#define _CRT_SECURE_NO_WARNINGS 1

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

/*------------------------------------------ VECTOR CLASS ------------------------------------------*/

// Implement path tracer TD1 

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
    Vector& operator+=(const Vector& V) {
        data[0] += V[0]; 
        data[1] += V[1]; 
        data[2] += V[2]; 
        return *this ;
    }
    
    double norm2() {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() {
		return sqrt(norm2());
	}
    void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
    const double& operator []( int i) const { return data[i] ;}
    double& operator []( int i) { return data[i]; }
    double data[3] ;
};

Vector operator+(const Vector& V, const Vector& U) {
    return Vector(V[0] + U[0] , V[1] + U[1] , V[2] + U[2]) ;
}
Vector div(Vector b, Vector a)
{
    return Vector(b[0] / a[0], b[1] / a[1], b[2] / a[2]);
}
Vector operator-(const Vector& V, const Vector& U){
    return Vector(V[0] - U[0] , V[1] - U[1] , V[2] - U[2]) ;
}
Vector operator-(const Vector& U){
    return Vector(- U[0] ,- U[1] ,- U[2]) ;
}
Vector operator*(const double alpha, const Vector& V) {
	return Vector(alpha*V[0], alpha*V[1], alpha*V[2]);
}
Vector operator*(const Vector& V, const double alpha) {
	return Vector(alpha*V[0], alpha*V[1], alpha*V[2]);
}
Vector operator/(const Vector& V, const double alpha) {
	return V*(1/alpha);
}

Vector cross(const Vector& V, const Vector &U) {
    return Vector(V[1] * U[2] - V[2] * U[1], V[2] * U[0] - V[0] * U[2], V[0] * U[1] - V[1] * U[0]);
}

double dot( const Vector& V, const Vector& U) { 
    return V[0] * U[0] + V[1] * U[1] + V[2] * U[2];
}

Vector operator*(const Vector& V, const Vector& U) {
	return Vector(V[0]*U[0], V[1]*U[1], V[2]*U[2]);       
}

Vector pow( Vector& V,  double n) {
    V[0]= std::pow(V[0], n);
    V[1]= std::pow(V[1], n);
    V[2]= std::pow(V[2], n);
	return Vector(V[0], V[1], V[2]);       
}

// Vector normalize(const Vector& V){
//     Vector U = V / V.norm();
//     return U ;
// }

