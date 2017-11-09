/*
 * write_Pontos.cpp
 *
 *  Created on: 25 de out de 2017
 *      Author: joaonizer
 *
 *      Describes the applied field as a 3D vctor
 */

#include "functions_sLLGS.h"
#include <iostream>

Eigen::RowVector3d compute_Happ(
        int step,
        int n_of_particles
        ){

    int size_of_s = 3; // Number of specificaton lines used in 's'
    int sim_step_index;  // This variable acumulates the steps already used to describe the applied field
    // its maximum value is going to be equals to N -> total number of simulation steps
    double slope;

    RowVector7d s[n_of_particles][size_of_s];
    Eigen::RowVector3d h_app;


    //=========================================//
    //======== Experiment Description =========//
    //=========================================//
    //
    //  s[particle][line] << xi, yi, zi, xf, yf, zf, N
    // #i - initial value of h_app on direction #
    // #f - final value of h_app on direction #
    // N - number os sim steps the variation occurs
    //

    s[0][0] << 0, 0, 0, 1, 1, 1, 10;
    s[0][1] << 1, 1, 1, 1, 1, 1, 10;
    s[0][2] << 1, 1, 1, 0, 0, 0, 10;

    s[1][0] << 0, 0, 0, -1, -1, -1, 10;
    s[1][1] << -1, -1, -1, -1, -1, -1, 10;
    s[1][2] << -1, -1, -1, 0, 0, 0, 10;

    //===== End of Experiment Description =====//



    for (int j=0; j<n_of_particles;j++){ // iterate on number of particles

        sim_step_index = 0; // Start index at 0

        for (int l=0; l<size_of_s;l++){ // iterate on lines of 's'

            for (int d = Coordinate::x ; d<=Coordinate::z; d++){ // iterate over the dimension (x,y,z)

                slope = (s[j][l](d+3)-s[j][l](d))/s[j][l](6); // Calculates the slope (x_f -x_i)/n;

                for (int n=0; s[j][l](d); n++){ // iterate from 0 to N

                    //h_app[sim_step_index+n][j](d) = s[j][l](d) +n*slope; // Caculates h_app = h_app_initial + n*slope

                } //end iterate on N
            } // end iterate on dimensions
            sim_step_index += s[j][l](6); // Increment index base on the number the line 'l' lasts
        } // end iterate on lines of 's'
    } //end iterate on particles

    //return h_app; // Return applied field
}
