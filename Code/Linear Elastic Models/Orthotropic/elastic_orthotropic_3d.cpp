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
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/elastic_orthotropic_3d.h"
#include "includes/checks.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

ElasticOrthotropic3D::ElasticOrthotropic3D()
    : ConstitutiveLaw()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

ElasticOrthotropic3D::ElasticOrthotropic3D(const ElasticOrthotropic3D& rOther)
    : ConstitutiveLaw(rOther)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer ElasticOrthotropic3D::Clone() const
{
    return Kratos::make_shared<ElasticOrthotropic3D>(*this);
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

ElasticOrthotropic3D::~ElasticOrthotropic3D()
{
};

/***********************************************************************************/
/***********************************************************************************/

void  ElasticOrthotropic3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;
    // b.- Get Values to compute the constitutive law:
    Flags & r_constitutive_law_options = rValues.GetOptions();

    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if( r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        CalculateCauchyGreenStrain( rValues, r_strain_vector);
    }

    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS )) {
        Vector& r_stress_vector = rValues.GetStressVector();
        CalculatePK2Stress( r_strain_vector, r_stress_vector, rValues);
    }

    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR )) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix( r_constitutive_matrix, rValues);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

// NOTE: Since we are in the hypothesis of small strains we can use the same function for everything

void ElasticOrthotropic3D::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

double& ElasticOrthotropic3D::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue)
{
    Vector& r_strain_vector = rParameterValues.GetStrainVector();
    Vector& r_stress_vector = rParameterValues.GetStressVector();

    if (rThisVariable == STRAIN_ENERGY) {
        this->CalculateCauchyGreenStrain(rParameterValues, r_strain_vector);
        this->CalculatePK2Stress( r_strain_vector, r_stress_vector, rParameterValues);

        rValue = 0.5 * inner_prod( r_strain_vector, r_stress_vector); // Strain energy = 0.5*E:C:E
    }

    return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

Vector& ElasticOrthotropic3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == STRAIN ||
        rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {
        this->CalculateCauchyGreenStrain( rParameterValues, rValue);
    } else if (rThisVariable == STRESSES ||
        rThisVariable == CAUCHY_STRESS_VECTOR ||
        rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
        rThisVariable == PK2_STRESS_VECTOR) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        // We compute the stress
        ElasticOrthotropic3D::CalculateMaterialResponseCauchy(rParameterValues);
        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    }

    return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& ElasticOrthotropic3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == CONSTITUTIVE_MATRIX ||
        rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        this->CalculateElasticMatrix(rValue, rParameterValues);
    }

    return( rValue );
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void ElasticOrthotropic3D::GetLawFeatures(Features& rFeatures)
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

int ElasticOrthotropic3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS_X);
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS_X] < 0.0) << "YOUNG_MODULUS_X is negative." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS_Y);
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS_Y] < 0.0) << "YOUNG_MODULUS_Y is negative." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS_Z);
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS_Z] < 0.0) << "YOUNG_MODULUS_Z is negative." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(SHEAR_MODULUS_XY);
    KRATOS_ERROR_IF(rMaterialProperties[SHEAR_MODULUS_XY] < 0.0) << "SHEAR_MODULUS_XY is negative." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(SHEAR_MODULUS_YZ);
    KRATOS_ERROR_IF(rMaterialProperties[SHEAR_MODULUS_YZ] < 0.0) << "SHEAR_MODULUS_YZ is negative." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(SHEAR_MODULUS_XZ);
    KRATOS_ERROR_IF(rMaterialProperties[SHEAR_MODULUS_XZ] < 0.0) << "SHEAR_MODULUS_XZ is negative." << std::endl;

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
	
	
	KRATOS_CHECK_VARIABLE_KEY(POISSON_RATIO_YZ);
	const double nu_3 = rMaterialProperties[POISSON_RATIO_YZ];
	KRATOS_ERROR_IF((nu_upper_bound - nu_3) < tolerance) << "POISSON_RATIO_YZ is above the upper bound 0.5." << std::endl;
    KRATOS_ERROR_IF((nu_3 - nu_lower_bound) < tolerance) << "POISSON_RATIO_YZ is below the lower bound -1.0." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0) << "DENSITY is negative." << std::endl;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::CheckClearElasticMatrix(Matrix& rConstitutiveMatrix)
{
    const SizeType size_system = this->GetStrainSize();
    if (rConstitutiveMatrix.size1() != size_system || rConstitutiveMatrix.size2() != size_system)
        rConstitutiveMatrix.resize(size_system, size_system, false);
    rConstitutiveMatrix.clear();
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::CalculateElasticMatrix(
    Matrix& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E1 = r_material_properties[YOUNG_MODULUS_X];
    const double E2 = r_material_properties[YOUNG_MODULUS_Y];
	const double E3 = r_material_properties[YOUNG_MODULUS_Z];
	
    const double v12 = r_material_properties[POISSON_RATIO_XY];
	const double v23 = r_material_properties[POISSON_RATIO_YZ];
	const double v13 = r_material_properties[POISSON_RATIO_XZ];
	
	const double G12 = r_material_properties[SHEAR_MODULUS_XY];
	const double G13 = r_material_properties[SHEAR_MODULUS_XZ];
	const double G23 = r_material_properties[SHEAR_MODULUS_YZ];
	
	this->CheckClearElasticMatrix(rConstitutiveMatrix);
	
	const double v21 =(E2*v12)/E1;
	const double v31 =(E3*v13)/E1;	
	const double v32 =(E3*v23)/E2;	

    const double Upsilon = (1.0)/(1.0-v12*v21-v23*v32-v13*v31-2.0*v21*v32*v13);

	rConstitutiveMatrix(0, 0) =  E1*Upsilon*(1.0-v23*v32);
	rConstitutiveMatrix(0, 1) =  E1*Upsilon*(v21+v31*v23);
	rConstitutiveMatrix(0, 2) =  E1*Upsilon*(v31+v21*v32);
	rConstitutiveMatrix(1, 0) =  E1*Upsilon*(v21+v31*v23);
	rConstitutiveMatrix(1, 1) =  E2*Upsilon*(1.0-v13*v31);
	rConstitutiveMatrix(1, 2) =  E2*Upsilon*(v32+v12*v31);
	rConstitutiveMatrix(2, 0) =  E1*Upsilon*(v31+v21*v32);
	rConstitutiveMatrix(2, 1) =  E2*Upsilon*(v32+v12*v31);
	rConstitutiveMatrix(2, 2) =  E3*Upsilon*(1.0-v12*v21);
	rConstitutiveMatrix(3, 3) =  G12;
	rConstitutiveMatrix(4, 4) =  G13;
	rConstitutiveMatrix(5, 5) =  G23;


}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::CalculatePK2Stress(
    const Vector& rStrainVector,
    Vector& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E1 = r_material_properties[YOUNG_MODULUS_X];
    const double E2 = r_material_properties[YOUNG_MODULUS_Y];
	const double E3 = r_material_properties[YOUNG_MODULUS_Z];
	
    const double v12 = r_material_properties[POISSON_RATIO_XY];
	const double v23 = r_material_properties[POISSON_RATIO_YZ];
	const double v13 = r_material_properties[POISSON_RATIO_XZ];
	
	const double G12 = r_material_properties[SHEAR_MODULUS_XY];
	const double G13 = r_material_properties[SHEAR_MODULUS_XZ];
	const double G23 = r_material_properties[SHEAR_MODULUS_YZ];
	

	const double v21 =(E2*v12)/E1;
	const double v31 =(E3*v13)/E1;	
	const double v32 =(E3*v23)/E2;	

    const double Upsilon = (1.0)/(1.0-v12*v21-v23*v32-v13*v31-2.0*v21*v32*v13);

	const double c11 =  E1*Upsilon*(1.0-v23*v32);
	const double c12 =  E1*Upsilon*(v21+v31*v23);
	const double c13 =  E1*Upsilon*(v31+v21*v32);
	const double c21 =  E1*Upsilon*(v21+v31*v23);
	const double c22 =  E2*Upsilon*(1.0-v13*v31);
	const double c23 =  E2*Upsilon*(v32+v12*v31);
	const double c31 =  E1*Upsilon*(v31+v21*v32);
	const double c32 =  E2*Upsilon*(v32+v12*v31);
	const double c33 =  E3*Upsilon*(1.0-v12*v21);
	const double c44 =  G12;
	const double c55 =  G13;
	const double c66 =  G23;
	

    rStressVector[0] =  c11 * rStrainVector[0] + c12 * rStrainVector[1] + c13 * rStrainVector[2];
    rStressVector[1] =  c21 * rStrainVector[0] + c22 * rStrainVector[1] + c23 * rStrainVector[2];
    rStressVector[2] =  c31 * rStrainVector[0] + c32 *  rStrainVector[1] + c33 * rStrainVector[2];
    rStressVector[3] = c44 * rStrainVector[3];
    rStressVector[4] = c55 * rStrainVector[4];
    rStressVector[5] = c66 * rStrainVector[5];
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticOrthotropic3D::CalculateCauchyGreenStrain(
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
