/* ************************************************************************************
 * Technical University of Munich 
 * Chair of Computational Mechanics
 * Progressive Damage Analysis of Laminated Composites (PDALAC) C++ Edition 
 * *************************************************************************************
 * Description: This tool is meant for testing the implementation of failure and damage models for composite laminates
 * Base constitutive law: Orthotropic 3D Linear Elasticity Model.
 * Failure Theories: Max Stress/Strain , Tsai-Wu, Hashin, Hashin-Rotem, Hoffman , Puck
 * Damage Theories: Ply Discounting Approach, Continum Damage Mechanics Approach
 * More Details: https://gitlab.lrz.de/Ammar.Kh/softwarelab19-group11-umat
 *  *************************************************************************************
 *  Author: Ammar Khallouf {Sept 2020}
 *  **************************************************************************************/
#include <stdio.h>
#include "ortho3d.h"

int main(int argc, char **argv)
{

    //****************************** Input Panel ******************************//

    //1. Declaring and Reading Orthotropic Material Elastic Properties Array ( E1,E2,E3,G12,G13,,G23,NU12,NU13,NU23)

    double input_mat_props[9]{164, 8.98, 8.98, 5.02, 5.02, 3.0, 0.32, 0.32, 0.496};

    //2. Declaring and Reading Orthotropic Material Strength Array ( Xt, Xc, Yt, Yc, Zt, Zc, , S12, S13, S23)

    double input_strength_props[9]{2.9, 1.68, 0.1, 0.247, 0.1, 0.247, 0.08, 0.08, 0.08};

    //3. Declaring and Reading Orthotropic Material Degredation Factors (Ply-Discount Approach)

    /*  beta_ft      degradation factor (fiber tension)           
     beta_fc      degradation factor (fiber compression)     
     beta_mt    degradation factor (matrix tension)          
     beta_mc    degradation factor (matrix compression)    
     beta_s       degradation factor (shear) */

    double input_dgrd_props[5]{0.5, 0.5, 0.5, 0.5, 0.5};

    //4. Declaring and Reading Strain Increments vector

    double input_strain_inc[6]{0.01, 0, 0, 0, 0, 0};

    //5. Declaring and Reading Number of Solution Steps

    size_t num_steps{6};

    //****************************** Analysis Panel ******************************//

    // 1.Creating (gauss_point) as an object of (ortho3d) type

    ortho3d gauss_point;

    // 2. Form 3D Orthotropic Elasticity Matrix

    gauss_point.clac_elastic_matrix(input_mat_props);

    // 3. Start Analysis Phase

    cout << "##### Starting Incremental Analysis #####" << endl;
    cout << "============================" << endl;

    for (size_t j = 0; j <= num_steps; j++)
    {
        cout << "\nStep No: " << j + 1 << endl;

        // Increment Strain Vector

        gauss_point.calc_strain_vector(input_strain_inc);

        //Calculate Stress Vector

        gauss_point.calc_stress_vector();

        //Calculate Failure indicies (Load/Strength ratios)

        gauss_point.calc_failure_vector(input_strength_props);

        //Calculate Damage, Modify Elasticity Matrix and Update Stresses accordingly->(if failure is present)

        gauss_point.calc_damage_vector(input_dgrd_props, input_strength_props);
    }

    return 0;
}
