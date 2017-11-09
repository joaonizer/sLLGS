/*
 * dm.cpp
 *
 *  Created on: 25 de out de 2017
 *      Author: joaonizer
 */


#include "functions_sLLGS.h"

Eigen::RowVector3d dm(
		Eigen::RowVector3d m,		// Magnetization of the particle
		Eigen::RowVector3d h_eff,	// Effective field of the particle
		Eigen::RowVector3d i_s,		// Spin current applied to the particle
		float alpha,				// Gilbert damping constant
		float alpha_l){			// Normalizing factor!!

	Eigen::RowVector3d m_is;		//	Cross Product m x i_s
	Eigen::RowVector3d mm_is;		//	Cross Product m x m_is
	Eigen::RowVector3d m_heff;		//	Cross Product m x h_eff
	Eigen::RowVector3d mm_heff;		//	Cross Product m x m_heff
	Eigen::RowVector3d dm;			//	Return variable dm

	m_is = m.cross(i_s);

    mm_is = m.cross(m_is);

	m_heff = m.cross(h_eff);

	mm_heff = m.cross(m_heff);

	return dm = -alpha_l*(m_heff+mm_is+alpha*(mm_heff-m_is));

}

