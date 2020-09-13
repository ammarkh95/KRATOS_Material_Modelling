#ifndef ORTHO3D
#define ORTHO3D

#include <iostream>
#include <vector>
#include <Eigen/Eigen/Dense>
#include <Eigen/Eigen/Sparse>

using namespace std;

using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6D;
typedef Matrix<double, 6, 1> Vector6D;

class ortho3d
{
private:

//ATTRIBUTES

Matrix6D Constitutive_Tensor=Matrix6D::Zero();
Vector6D Stress_Vector=Vector6D::Zero();
Vector6D Strain_Vector=Vector6D::Zero();
Vector6D Failure_indicies=Vector6D::Zero();
Vector6D Damage_Coeffs=Vector6D::Ones();


public:

//default constructor 

ortho3d();

//destructor

~ortho3d();

//methods decleration

void clac_elastic_matrix (double m_props[9]);

void calc_stress_vector ();

void calc_strain_vector (double strain_inc[6] );

void calc_failure_vector (double strength[9]);

void calc_damage_vector (double dgrd[5], double strength[9]);


};

#endif // ORTHO3D
