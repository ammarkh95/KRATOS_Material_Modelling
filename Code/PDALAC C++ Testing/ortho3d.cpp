#include "ortho3d.h"

//constructor

ortho3d::ortho3d()
{
    //    cout << "constructor is called" << endl;
}
//**************************************************************************************//

//destructor
ortho3d::~ortho3d()
{

    //    cout << "destructor is called" << endl;
}

//**************************************************************************************//

//Form Material Elasticity Matrix

void ortho3d::clac_elastic_matrix(double m_props[9])
{

    //Read Mechanics Properites of Orthotropic Material (9 Constants)

    const double E11 = m_props[0];
    const double E22 = m_props[1];
    const double E33 = m_props[2];
    const double G12 = m_props[3];
    const double G13 = m_props[4];
    const double G23 = m_props[5];
    const double NU12 = m_props[6];
    const double NU13 = m_props[7];
    const double NU23 = m_props[8];

    //Calculate the conjugate poisson ratios from symmetry relations

    const double NU21 = (E22 * NU12) / E11;
    const double NU31 = (E33 * NU13) / E11;
    const double NU32 = (E33 * NU23) / E22;

    const double Upsilon = (1.0) / (1.0 - NU12 * NU21 - NU23 * NU32 - NU13 * NU31 - 2.0 * NU21 * NU32 * NU13);

    //Fill the coefficents of the Elasticity Matrix as per ABAQUS Convention

    Constitutive_Tensor(0, 0) = E11 * Upsilon * (1.0 - NU23 * NU32);
    Constitutive_Tensor(0, 1) = E11 * Upsilon * (NU21 + NU31 * NU23);
    Constitutive_Tensor(0, 2) = E11 * Upsilon * (NU31 + NU21 * NU32);
    Constitutive_Tensor(1, 0) = E11 * Upsilon * (NU21 + NU31 * NU23);
    Constitutive_Tensor(1, 1) = E22 * Upsilon * (1.0 - NU13 * NU31);
    Constitutive_Tensor(1, 2) = E22 * Upsilon * (NU32 + NU12 * NU31);
    Constitutive_Tensor(2, 0) = E11 * Upsilon * (NU31 + NU21 * NU32);
    Constitutive_Tensor(2, 1) = E22 * Upsilon * (NU32 + NU12 * NU31);
    Constitutive_Tensor(2, 2) = E33 * Upsilon * (1.0 - NU12 * NU21);
    Constitutive_Tensor(3, 3) = G12;
    Constitutive_Tensor(4, 4) = G13;
    Constitutive_Tensor(5, 5) = G23;

    cout << "Elasticity Matrix Successfully Formed\n";
    cout << "=======================\n";
}

//**************************************************************************************//

//Calculate Strain Vector by acummulating strain increments

void ortho3d::calc_strain_vector(double strain_inc[6])
{
    for (size_t i = 0; i < 6; i++)
    {

        Strain_Vector(i) += strain_inc[i];
    }

    //Output

    IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
    cout << "\nStrain Vector: \n";
    cout << Strain_Vector.format(CommaInitFmt) << endl;
    cout << "=======================\n";
}

//**************************************************************************************//

//Calculate Stress Vector

void ortho3d::calc_stress_vector()
{

    Stress_Vector = Constitutive_Tensor * Strain_Vector;

    //Output
    IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
    cout << "Stress Vector: \n";
    cout << Stress_Vector.format(CommaInitFmt) << endl;
    cout << "=======================\n";
}

//**************************************************************************************//

//Evaluate Failure Indicies

void ortho3d::calc_failure_vector(double strength[9])
{

    //Fiber Failure Tension/Compression

    if (Stress_Vector(0) >= 0)
    {
        Failure_indicies(0) = Stress_Vector(0) / strength[0];
    }
    else
    {
        Failure_indicies(0) = Stress_Vector(0) / strength[1];
    }

    //Matrix Failure Tension/Compression

    if (Stress_Vector(1) >= 0)
    {
        Failure_indicies(1) = Stress_Vector(1) / strength[2];
    }
    else
    {
        Failure_indicies(1) = Stress_Vector(1) / strength[3];
    }

    //Interlaminar Failure Tension/Compression

    if (Stress_Vector(2) >= 0)
    {
        Failure_indicies(2) = Stress_Vector(2) / strength[4];
    }

    else
    {
        Failure_indicies(2) = Stress_Vector(2) / strength[5];
    }

    //Shear Failure

    Failure_indicies(3) = Stress_Vector(3) / strength[6];

    Failure_indicies(4) = Stress_Vector(4) / strength[7];

    Failure_indicies(5) = Stress_Vector(5) / strength[8];

    //Output
    IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
    cout << "Failure Indicies: \n";
    cout << Failure_indicies.format(CommaInitFmt) << endl;
    cout << "=======================\n";
}

//**************************************************************************************//

//Evaluate Damage Coefficents, Modify Elasticity Matrix, Update Stresses ->(if failure is present)

void ortho3d::calc_damage_vector(double dgrd[5], double strength[9])
{

    // While failure is detected apply damage and update stresses

    while (Failure_indicies.maxCoeff() > 1.0)
    {

        //Fiber Damage Tension/Compression

        if (Stress_Vector(0) >= 0 and Failure_indicies(0) > 1)
        {
            Damage_Coeffs(0) = abs(0.9999 - dgrd[0]);
        }
        else if (Stress_Vector(0) < 0 and Failure_indicies(0) > 1)
        {
            Damage_Coeffs(0) = abs(0.9999 - dgrd[1]);
        }

        //Matrix Damage Tension/Compression

        if (Stress_Vector(1) >= 0 and Failure_indicies(1) > 1)
        {
            Damage_Coeffs(1) = abs(0.9999 - dgrd[2]);
        }

        else if (Stress_Vector(1) < 0 and Failure_indicies(1) > 1)
        {
            Damage_Coeffs(1) = abs(0.9999 - dgrd[3]);
        }

        //Interlaminar Damage Tension/Compression

        if (Stress_Vector(2) >= 0 and Failure_indicies(2) > 1)
        {
            Damage_Coeffs(2) = abs(0.9999 - dgrd[2]);
        }

        else if (Stress_Vector(2) < 0 and Failure_indicies(2) > 1)
        {
            Damage_Coeffs(2) = abs(0.9999 - dgrd[3]);
        }

        //Shear Damage

        for (size_t j = 3; j <= 5; j++)
        {
            if (Failure_indicies(j) > 1)
            {
                Damage_Coeffs(j) = abs(0.9999 - dgrd[4]);
                cout << "ds is " << Damage_Coeffs(j) << endl;
            }
        }

        //Factorize the Elasticity Matrix Coefficents with Damage Coefficents

        Constitutive_Tensor(0, 0) *= Damage_Coeffs(0);
        Constitutive_Tensor(0, 1) *= Damage_Coeffs(0) * Damage_Coeffs(1);
        Constitutive_Tensor(0, 2) *= Damage_Coeffs(0) * Damage_Coeffs(2);
        Constitutive_Tensor(1, 0) *= Damage_Coeffs(0) * Damage_Coeffs(1);
        Constitutive_Tensor(1, 1) *= Damage_Coeffs(1);
        Constitutive_Tensor(1, 2) *= Damage_Coeffs(1) * Damage_Coeffs(2);
        Constitutive_Tensor(2, 0) *= Damage_Coeffs(0) * Damage_Coeffs(2);
        Constitutive_Tensor(2, 1) *= Damage_Coeffs(1) * Damage_Coeffs(2);
        Constitutive_Tensor(2, 2) *= Damage_Coeffs(2);
        Constitutive_Tensor(3, 3) *= Damage_Coeffs(3);
        Constitutive_Tensor(4, 4) *= Damage_Coeffs(4);
        Constitutive_Tensor(5, 5) *= Damage_Coeffs(5);

        //Output
        IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
        cout << " Damage Coefficents: \n";
        cout << Damage_Coeffs.format(CommaInitFmt) << endl;
        cout << "=======================\n";

        //Update Stress Vector using Modified Elasticity Matrix
        calc_stress_vector();
        //Update Failure indicies
        calc_failure_vector(strength);
    }
}