/*
 * compute_Heff.cpp
 *
 *  Created on: 24 de out de 2017
 *      Author: joaonizer
 *
 *      This function calculates the effective field;
 */

#include "functions_sLLGS.h"


Eigen::RowVector3d compute_Heff(
        Eigen::RowVector3d m,		// Magnetization of the particle
        Eigen::RowVector3d h_app,	// Applied field at that instant of time
        double hkms,				// Anisotropy field normalized by Ms
        Eigen::Matrix3d Nd,			// Demagnetizing Tensor
        Eigen::RowVector3d h_c,		// Coupling Field
        Eigen::RowVector3d h_th,	// Thermal field
        Eigen::RowVector3d n) {		// Easy axis of the particle.

    Eigen::RowVector3d h_eff;  		// Return variable of effective field

    double nm = m*n.adjoint();		// Dot product of n and m

    Eigen::RowVector3d h_ans = hkms * nm * n; 		// Anisotropy field

    Eigen::RowVector3d h_d = m * Nd.transpose(); 	// Demagnetizing field

    //h_eff = h_app + h_ans - h_d -h_c + h_th;		// Calculation of the effective field
    h_eff = h_app + h_ans - h_d -h_c + h_th;		// removed hc for debub
    return h_eff;
}

