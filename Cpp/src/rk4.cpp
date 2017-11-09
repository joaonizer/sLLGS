/*
 * rk4.cpp
 *
 *  Created on: 25 de out de 2017
 *      Author: joaonizer
 *
 *      Runge-Kuta method for solving LLGS
 */

#include "functions_sLLGS.h"


Eigen::RowVector3d rk4(
		Eigen::RowVector3d m,
		Eigen::RowVector3d h_eff,
		Eigen::RowVector3d i_s,
		double dt,				// Normalized time_step (non-dimensional)
		float alpha,			// Gilbert damping constant
		float alpha_l){			// Normalizing factor!!


	Eigen::RowVector3d k[4];
	Eigen::RowVector3d mm[3];
	Eigen::RowVector3d m_new;

	k[0]  = dm(m, h_eff, i_s, alpha, alpha_l);
	mm[0] = m+k[0]*dt/2;

	k[1]  = dm(mm[0], h_eff, i_s, alpha, alpha_l);
	mm[1] = m+k[1]*dt/2;

	k[2]  = dm(mm[1], h_eff, i_s, alpha, alpha_l);
	mm[2] = m+k[2]*dt;

	k[3]  = dm(mm[2], h_eff, i_s, alpha, alpha_l);

    m_new = m+(k[0] + 2*(k[1] + k[2]) + k[3])*dt/6;

	return m_new;

}



