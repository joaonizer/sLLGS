/*
 * write_Pontos.cpp
 *
 *  Created on: 25 de out de 2017
 *      Author: joaonizer
 */

#include "functions_sLLGS.h"
#include <iostream>


void write_Points(
        double *w,
        double *l,
        double dx,
        Eigen::RowVector3d *d_or,
        Eigen::RowVector4d cortes_y,
        int j,
        Eigen::RowVector4d *px,
        Eigen::RowVector4d *py
        ){

    // Compute points X
    double x1 = -w[j]/2;
    double x2 =  w[j]/2;

    // Computes points Y
    double y1 =  l[j]/2;
    double y2 = -l[j]/2;

    //std::cout << w[j] << std::endl;
    //std::cout << l[j] << std::endl;

    // Create vector px_j and py_j
    px[j] = {x1, x2, x2, x1};
//    py[j] = {y1, y1, y2, y2} - cortes_y; // For slanted shapes
    py[j] << y1, y1, y2, y2; // For slanted shapes
    py[j] -= cortes_y;

    // Computes distance of origin for particle j
    if (j==0)
        d_or[j] << 0,0,0;
    else {
        d_or[j] << (w[j-1]/2 + w[j]/2 + dx + d_or[j-1](Coordinate::x)) ,0 ,0;
    //std::cout << j-1 << " " << d_or[j-1] << " " << j << " " << d_or[j] << std::endl;
    }
    // adds w/2 of particle j-1
    // adds w/2 of particle j
    // adds dx - distance between particle j-1 and j
    // adds cumulative X shift on particle j-1

}
