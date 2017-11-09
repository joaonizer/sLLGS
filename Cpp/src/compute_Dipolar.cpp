/*
 * compute_Dipolar.cpp
 *
 *  Created on: 25 de out de 2017
 *      Author: joaonizer
 *      Reference:  dipolar3D3.f90              Author: Luiz G. Melo
 *                  dipolar3D3_modificado.f90   Editado: João Nizer
 */

#include "functions_sLLGS.h"
#include <iostream>


extern"C" {
void dipolar3D3(
        double *ppx_j,
        double *ppy_j,
        double *th_j,
        double *ddo_j,
        double *ppx_k,
        double *ppy_k,
        double *th_k,
        double *ddo_k,
        double *pc
        );
}

Eigen::Matrix3d compute_Dipolar(
        Eigen::RowVector4d px_j,    // X-positions of particle j
        Eigen::RowVector4d py_j,    // Y-positions of particle j
        double th_j,                 // Thickness of particle j
        Eigen::RowVector3d do_j,    // Center of particle j
        Eigen::RowVector4d px_k,    // X-positions of particle k
        Eigen::RowVector4d py_k,    // Y-positions of particle k
        double th_k,                 // Thickness of particle  k
        Eigen::RowVector3d do_k     // Center of particle k
        ){

    double ppx_j [4] = {px_j(0), px_j(1), px_j(2), px_j(3)};
    double ppy_j [4] = {py_j(0), py_j(1), py_j(2), py_j(3)};
    double ddo_j [3] = {do_j(0), do_j(1), do_j(2)};

    double ppx_k [4] = {px_k(0), px_k(1), px_k(2), px_k(3)};
    double ppy_k [4] = {py_k(0), py_k(1), py_k(2), py_k(3)};
    double ddo_k [3] = {do_k(0), do_k(1), do_k(2)};

    double pc [7] = { 0 };
    Eigen::Matrix3d ppc;

    dipolar3D3(ppx_j,ppy_j,&th_j,ddo_j,ppx_k,ppy_k,&th_k,ddo_k,pc);

    // Print out the components
    //for (double p:pc)
    //std::cout << p << std::endl;

    // Defines the Tensor based 'pc'
    ppc(0,0) = pc[0] + pc[1];
    ppc(0,1) = pc[4];
    ppc(0,2) = pc[5];

    ppc(1,0) = pc[4];
    ppc(1,1) = pc[0] + pc[2];
    ppc(1,2) = pc[6];

    ppc(2,0) = pc[5];
    ppc(2,1) = pc[6];
    ppc(2,2) = pc[0] + pc[3];

    return ppc; //Return the Tensor

}

