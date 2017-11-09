/*
 * compute_Demag.cpp
 *
 *  Created on: 25 de out de 2017
 *      Author: joaonizer
 *      Reference:  demag3D3.f90              Author: Luiz G. Melo
 *                  demag3D3_modificado.f90   Editado: João Nizer
 */

#include "functions_sLLGS.h"
#include <iostream>


extern"C" {
void demag3D3(double *px, double *py, double *th, double *d);
}

Eigen::Matrix3d compute_Demag(
        Eigen::RowVector4d px,
        Eigen::RowVector4d py,
        double th
        ){

	double ppx [4] = {px(0), px(1), px(2), px(3)};
	double ppy [4] = {py(0), py(1), py(2), py(3)};
	double pd [9] = { 0 };
	Eigen::Matrix3d ppd;

    demag3D3(ppx,ppy,&th,pd);

	ppd(0,0) = pd[0];
	ppd(0,1) = pd[1];
	ppd(0,2) = pd[2];
	ppd(1,0) = pd[3];
	ppd(1,1) = pd[4];
	ppd(1,2) = pd[5];
	ppd(2,0) = pd[6];
	ppd(2,1) = pd[7];
	ppd(2,2) = pd[8];

	return ppd;

}

