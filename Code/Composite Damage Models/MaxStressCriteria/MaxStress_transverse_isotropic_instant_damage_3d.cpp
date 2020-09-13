// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//
//  Contributor :    Ammar Khallouf 
//
// System includes

#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/MaxStress_transverse_isotropic_instant_damage_3d.h"
#include "includes/checks.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

MaxStressTransverseIsotropicInstantDamage::MaxStressTransverseIsotropicInstantDamage()
    : ConstitutiveLaw()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

MaxStressTransverseIsotropicInstantDamage::MaxStressTransverseIsotropicInstantDamage(const MaxStressTransverseIsotropicInstantDamage& rOther)
    : ConstitutiveLaw(rOther)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer MaxStressTransverseIsotropicInstantDamage::Clone() const
{
    return Kratos::make_shared<MaxStressTransverseIsotropicInstantDamage>(*this);
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

MaxStressTransverseIsotropicInstantDamage::~MaxStressTransverseIsotropicInstantDamage()
{
};

/***********************************************************************************/
/***********************************************************************************/

void  MaxStressTransverseIsotropicInstantDamage::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;
    // b.- Get Values to compute the constitutive law:
    Flags & r_constitutive_law_options = rValues.GetOptions();

    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if( r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        CalculateCauchyGreenStrain( rValues, r_strain_vector);
    }

    //Elastic Matrix Computed Here
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR )) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix( r_constitutive_matrix,rValues);
    }

    // We compute the stress
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS )) {

        //Elastic Matrix
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix( r_constitutive_matrix,rValues);

        //Calculate Predictive Stresses
        Vector& r_stress_vector = rValues.GetStressVector();
        CalculatePK2Stress( r_strain_vector, r_stress_vector,r_constitutive_matrix);

        //Check Failure as per Max Stress Criteria
        CheckFailureMaxStress(r_strain_vector,r_stress_vector,rValues);

        //Modify Constitutive Tensor based on detected Failure Mode
        CalculateSecantMatrix(r_constitutive_matrix,r_stress_vector, r_strain_vector,rValues);

        //Update Stresses
        CalculatePK2Stress( r_strain_vector, r_stress_vector,r_constitutive_matrix);        

    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

// NOTE: Since we are in the hypothesis of small strains we can use the same function for everything

void MaxStressTransverseIsotropicInstantDamage::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ANISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = 6;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}

/***********************************************************************************/
/***********************************************************************************/

int MaxStressTransverseIsotropicInstantDamage::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS_X);
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS_X] < 0.0) << "YOUNG_MODULUS_X is negative." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS_Y);
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS_Y] < 0.0) << "YOUNG_MODULUS_Y is negative." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(POISSON_RATIO_XY);
    const double tolerance = 1.0e-12;
    const double nu_upper_bound = 0.5;
    const double nu_lower_bound = -1.0;
    const double nu_1 = rMaterialProperties[POISSON_RATIO_XY];
    KRATOS_ERROR_IF((nu_upper_bound - nu_1) < tolerance) << "POISSON_RATIO_XY is above the upper bound 0.5." << std::endl;
    KRATOS_ERROR_IF((nu_1 - nu_lower_bound) < tolerance) << "POISSON_RATIO_XY is below the lower bound -1.0." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(POISSON_RATIO_XZ);
    const double nu_2 = rMaterialProperties[POISSON_RATIO_XZ];
    KRATOS_ERROR_IF((nu_upper_bound - nu_2) < tolerance) << "POISSON_RATIO_XZ is above the upper bound 0.5." << std::endl;
    KRATOS_ERROR_IF((nu_2 - nu_lower_bound) < tolerance) << "POISSON_RATIO_XZ is below the lower bound -1.0." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0) << "DENSITY is negative." << std::endl;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::CheckClearElasticMatrix(Matrix& rConstitutiveMatrix)
{
    const SizeType size_system = this->GetStrainSize();
    if (rConstitutiveMatrix.size1() != size_system || rConstitutiveMatrix.size2() != size_system)
        rConstitutiveMatrix.resize(size_system, size_system, false);
    rConstitutiveMatrix.clear();
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::CalculateElasticMatrix(
    Matrix& rConstitutiveMatrix,ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E_1 = r_material_properties[YOUNG_MODULUS_X];
    const double E_2 = r_material_properties[YOUNG_MODULUS_Y];
    const double NU_12 = r_material_properties[POISSON_RATIO_XY];
    const double NU_23 = r_material_properties[POISSON_RATIO_YZ];
    const double G_12   = r_material_properties[SHEAR_MODULUS_XY];

    this->CheckClearElasticMatrix(rConstitutiveMatrix);

    const double NU_21 = (E_2*NU_12)/E_1;
    const double G_23 = (E_2)/(2.0*(1.0+NU_23));
    const double Upsilon = 1.0/(1.0-NU_23*NU_23-2*NU_12*NU_21-2*NU_23*NU_21*NU_12);


    rConstitutiveMatrix( 0, 0 ) = E_1*Upsilon*(1.0-NU_23*NU_23);
    rConstitutiveMatrix( 0, 1 ) = E_1*Upsilon*(NU_21+ NU_21*NU_23);
    rConstitutiveMatrix( 0, 2 ) = E_1*Upsilon*(NU_21+ NU_21*NU_23);
    rConstitutiveMatrix( 1, 0 ) = E_1*Upsilon*(NU_21+ NU_21*NU_23);
    rConstitutiveMatrix( 1, 1 ) = E_2*Upsilon*(1.0- NU_12*NU_21);
    rConstitutiveMatrix( 1, 2 ) = E_2*Upsilon*(NU_23+ NU_12*NU_21);
    rConstitutiveMatrix( 2, 0 ) = E_1*Upsilon*(NU_21+ NU_21*NU_23);
    rConstitutiveMatrix( 2, 1 ) = E_2*Upsilon*(NU_23+ NU_12*NU_21);
    rConstitutiveMatrix( 2, 2 ) = E_2*Upsilon*(1.0- NU_12*NU_21);
    rConstitutiveMatrix( 3, 3 ) = G_12;
    rConstitutiveMatrix( 4, 4 ) = G_12;
    rConstitutiveMatrix( 5, 5 ) = G_23;

    //std::cout<<"Initial Elastic Matrix Formed"<<std::endl;


// Update Constiutive Tensor for damage from previous steps (if present)

    rConstitutiveMatrix( 0, 0 ) *= Damage_Coeffs[0];
    rConstitutiveMatrix( 0, 1 ) *= Damage_Coeffs[0]*Damage_Coeffs[1];
    rConstitutiveMatrix( 0, 2 ) *= Damage_Coeffs[0]*Damage_Coeffs[2];
    rConstitutiveMatrix( 1, 0 ) *= Damage_Coeffs[0]*Damage_Coeffs[1];
    rConstitutiveMatrix( 1, 1 ) *= Damage_Coeffs[1];
    rConstitutiveMatrix( 1, 2 ) *= Damage_Coeffs[1]*Damage_Coeffs[2];
    rConstitutiveMatrix( 2, 0 ) *= Damage_Coeffs[0]*Damage_Coeffs[2];
    rConstitutiveMatrix( 2, 1 ) *= Damage_Coeffs[1]*Damage_Coeffs[2];
    rConstitutiveMatrix( 2, 2 ) *= Damage_Coeffs[2];
    rConstitutiveMatrix( 3, 3 ) *= Damage_Coeffs[3];
    rConstitutiveMatrix( 4, 4 ) *= Damage_Coeffs[4];
    rConstitutiveMatrix( 5, 5 ) *= Damage_Coeffs[5];


}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::CalculatePK2Stress(
    const Vector& rStrainVector,Vector& rStressVector,Matrix& rConstitutiveMatrix)
{

    //Calculating Stress Tensor

    rStressVector[0] = rConstitutiveMatrix( 0, 0 )* rStrainVector[0]
     + rConstitutiveMatrix( 0, 1 )*rStrainVector[1] 
     +  rConstitutiveMatrix( 0, 2 )*rStrainVector[2];

    rStressVector[1] =  rConstitutiveMatrix( 1, 0 )* rStrainVector[0]
     + rConstitutiveMatrix( 1, 1 )*rStrainVector[1] 
     +  rConstitutiveMatrix( 1, 2 )*rStrainVector[2];

    rStressVector[2] =  rConstitutiveMatrix( 2, 0 )* rStrainVector[0]
      +  rConstitutiveMatrix( 2, 1 )* rStrainVector[1]
      + rConstitutiveMatrix( 2, 2 )*rStrainVector[2];

    rStressVector[3] =  rConstitutiveMatrix( 3, 3 )*rStrainVector[3];

    rStressVector[4] =  rConstitutiveMatrix( 4, 4 )*rStrainVector[4];

    rStressVector[5] = rConstitutiveMatrix( 5, 5 )* rStrainVector[5];


    //std::cout<<"CalculatePK2Stress: "<<rStressVector[0]<<std::endl;

  
}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::CheckFailureMaxStress(const Vector& rStrainVector,
    const Vector& rStressVector,ConstitutiveLaw::Parameters& rValues)
    
{

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    //Reading Laminate Strength Parameters
    const double Xt = r_material_properties[LAMINATE_TENSILE_STRENGTH_11];
    const double Xc = r_material_properties[LAMINATE_COMPRESSIVE_STRENGTH_11];
    const double Yt = r_material_properties[LAMINATE_TENSILE_STRENGTH_22];
    const double Yc = r_material_properties[LAMINATE_COMPRESSIVE_STRENGTH_22];
    const double Zt = r_material_properties[LAMINATE_INTERLAMINAR_TENSILE_STRENGTH_33];
    const double Zc = r_material_properties[LAMINATE_INTERLAMINAR_COMPRESSIVE_STRENGTH_33];
    const double S12 = r_material_properties[LAMINATE_SHEAR_STRENGTH_12];
    const double S13 = r_material_properties[LAMINATE_SHEAR_STRENGTH_13];
    const double S23 = r_material_properties[LAMINATE_SHEAR_STRENGTH_23];


    //Evaluate Failure indicies :

//1. Failure along fiber (11 axis) (Tension/Compression)

    if(rStressVector[0] >= 0.0){

       Failure_Indicies[0]= rStressVector[0]/Xt;
    }

    else {
       Failure_Indicies[0]= -rStressVector[0]/Xc;
    }


//2. Failure along matrix (22 axis) (Tension/Compression)

    if(rStressVector[1] >= 0.0){

       Failure_Indicies[1]= rStressVector[1]/Yt;
    }

    else {
       Failure_Indicies[1]= -rStressVector[1]/Yc;

    }

//3. Interlaminar failure (33 axis) (Tension/Compression)                   


    if(rStressVector[2] >= 0.0){

       Failure_Indicies[2]= rStressVector[2]/Zt;
    }

    else {
       Failure_Indicies[2]= -rStressVector[2]/Zc;
    }

//4. Shear Failure (1-2 plane)

    Failure_Indicies[3]= abs(rStressVector[3])/S12;


//5. Shear Failure (1-3 plane)

    Failure_Indicies[4]= abs(rStressVector[4])/S13;


//6. Shear Failure (2-3 plane)

    Failure_Indicies[5]= abs(rStressVector[5])/S23;


    //std::cout<<"CheckFailureMaxStress: "<<Failure_Indicies[0]<<std::endl;

}

/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::CalculateSecantMatrix(
        Matrix& rConstitutiveMatrix, Vector& rStressVector,const Vector& rStrainVector,
        ConstitutiveLaw::Parameters& rValues)
    
{

    //Damage Iterations Loop (Accessed only if failure is detected)

     while (*std::max_element(Failure_Indicies,Failure_Indicies+6)> 1.0)
    {


    const Properties& r_material_properties = rValues.GetMaterialProperties();

    //Reading Degredation Parameters

    /*/const double Xt = r_material_properties[LAMINATE_TENSILE_STRENGTH_11];
    const double Xc = r_material_properties[LAMINATE_COMPRESSIVE_STRENGTH_11];
    const double Yt = r_material_properties[LAMINATE_TENSILE_STRENGTH_22];
    const double Yc = r_material_properties[LAMINATE_COMPRESSIVE_STRENGTH_22];
    const double Zt = r_material_properties[LAMINATE_INTERLAMINAR_TENSILE_STRENGTH_33];
    const double Zc = r_material_properties[LAMINATE_INTERLAMINAR_COMPRESSIVE_STRENGTH_33];
    const double S12 = r_material_properties[LAMINATE_SHEAR_STRENGTH_12];
    const double S13 = r_material_properties[LAMINATE_SHEAR_STRENGTH_13];
    const double S23 = r_material_properties[LAMINATE_SHEAR_STRENGTH_23];*/
    const double beta_ft = r_material_properties[DAMAGE_FIBER_TENSION_11];
    const double beta_fc = r_material_properties[DAMAGE_FIBER_COMPRESSION_11];
    const double beta_mt = r_material_properties[DAMAGE_MATRIX_TENSION_22];
    const double beta_mc = r_material_properties[DAMAGE_MATRIX_COMPRESSION_22];
    const double beta_s = r_material_properties[DAMAGE_SHEAR];

    //Go through failure indicies and apply damage accordingly

//1. Damage along fiber (11 axis) (Tension/Compression)

    if(rStressVector[0] >= 0.0 and Failure_Indicies[0]>1.0){

       Damage_Coeffs[0]= fabs(0.99-beta_ft);
       //std::cout<<"CallingDamage: "<<Damage_Coeffs[0]<<std::endl;

    }

    else if (rStressVector[0] < 0.0 and Failure_Indicies[0]>1.0){
       Damage_Coeffs[0]= fabs(0.99-beta_fc);
    
    }

//2. Damage along matrix (22 axis) (Tension/Compression)

    if(rStressVector[1] >= 0.0 and Failure_Indicies[1]>1.0){

       Damage_Coeffs[1]= fabs(0.99-beta_mt);
    }

    else if (rStressVector[1] < 0.0 and Failure_Indicies[1]>1.0){
       Damage_Coeffs[1]= fabs(0.99-beta_mc);
    
    }

//3. Interlaminar Damage (33 axis) (Tension/Compression)                   


    if(rStressVector[2] >= 0.0 and Failure_Indicies[2]>1.0){

       Damage_Coeffs[2]= fabs(0.99-beta_mt);
    }

    else if (rStressVector[2] < 0.0 and Failure_Indicies[2]>1.0){
       Damage_Coeffs[2]= fabs(0.99-beta_mc);
    
    }

//4. Shear Damage (1-2 plane)

    if(Failure_Indicies[3]>1.0){

       Damage_Coeffs[3]= fabs(0.99-beta_s);
    }

//5. Shear Damage (1-3 plane)

    if(Failure_Indicies[4]>1.0){

       Damage_Coeffs[4]= fabs(0.99-beta_s);
    }

//6. Shear Damage (2-3 plane)

    if(Failure_Indicies[5]>1.0){

       Damage_Coeffs[5]= fabs(0.99-beta_s);
    }

    rConstitutiveMatrix( 0, 0 ) *= Damage_Coeffs[0];
    rConstitutiveMatrix( 0, 1 ) *= Damage_Coeffs[0]*Damage_Coeffs[1];
    rConstitutiveMatrix( 0, 2 ) *= Damage_Coeffs[0]*Damage_Coeffs[2];
    rConstitutiveMatrix( 1, 0 ) *= Damage_Coeffs[0]*Damage_Coeffs[1];
    rConstitutiveMatrix( 1, 1 ) *= Damage_Coeffs[1];
    rConstitutiveMatrix( 1, 2 ) *= Damage_Coeffs[1]*Damage_Coeffs[2];
    rConstitutiveMatrix( 2, 0 ) *= Damage_Coeffs[0]*Damage_Coeffs[2];
    rConstitutiveMatrix( 2, 1 ) *= Damage_Coeffs[1]*Damage_Coeffs[2];
    rConstitutiveMatrix( 2, 2 ) *= Damage_Coeffs[2];
    rConstitutiveMatrix( 3, 3 ) *= Damage_Coeffs[3];
    rConstitutiveMatrix( 4, 4 ) *= Damage_Coeffs[4];
    rConstitutiveMatrix( 5, 5 ) *= Damage_Coeffs[5];

    //Update Stresses
    this->CalculatePK2Stress(rStrainVector,rStressVector,rConstitutiveMatrix);

    //Update Failure Indicies (Check for further failure)
    this->CheckFailureMaxStress(rStrainVector,rStressVector,rValues);

    }

}


/***********************************************************************************/
/***********************************************************************************/

void MaxStressTransverseIsotropicInstantDamage::CalculateCauchyGreenStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector
    )
{
    const SizeType space_dimension = this->WorkingSpaceDimension();

    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();
    KRATOS_DEBUG_ERROR_IF(F.size1()!= space_dimension || F.size2() != space_dimension)
        << "expected size of F " << space_dimension << "x" << space_dimension << ", got " << F.size1() << "x" << F.size2() << std::endl;

    Matrix E_tensor = prod(trans(F),F);
    for(unsigned int i=0; i<space_dimension; ++i)
      E_tensor(i,i) -= 1.0;
    E_tensor *= 0.5;

    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

} // Namespace Kratos
