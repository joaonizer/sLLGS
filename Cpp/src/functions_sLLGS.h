/*
 * functions_sLLGS.h
 *
 *  Created on: 24 de out de 2017
 *      Author: joaonizer
 */

#ifndef FUNCTIONS_SLLGS_H_
#define FUNCTIONS_SLLGS_H_

#include <eigen3/Eigen/Dense>
#include <stdio.h>
#include <stdlib.h>
#include <string>

typedef enum {
    x = 0,
    y = 1,
    z = 2
} Coordinate;

// Function to Calculate the Effective Field //
Eigen::RowVector3d compute_Heff(
        Eigen::RowVector3d m,
        Eigen::RowVector3d h_app,
        double hkms,
        Eigen::Matrix3d Nd,
        Eigen::RowVector3d h_c,
        Eigen::RowVector3d h_th,
        Eigen::RowVector3d n
        );

// Function to calculate the variation of magnetization //
Eigen::RowVector3d dm(
        Eigen::RowVector3d m,
        Eigen::RowVector3d h_eff,
        Eigen::RowVector3d i_s,
        float alpha,
        float alpha_l
        );

// Range-Kuta Method for Solving s-LLGS
Eigen::RowVector3d rk4(
        Eigen::RowVector3d m,
        Eigen::RowVector3d h_eff,
        Eigen::RowVector3d i_s,
        double dt,
        float alpha,
        float alpha_l
        );

// Function to Calculate the Demagnetizing Tensor //
Eigen::Matrix3d compute_Demag(
        Eigen::RowVector4d px,      // X-positions of partcle
        Eigen::RowVector4d py,      // Y-positions of partcle
        double th                   // Thickness of particle
        );

// Function to Calculate the Caoupling Tensor
Eigen::Matrix3d compute_Dipolar(
        Eigen::RowVector4d px_j,    // X-positions of particle j
        Eigen::RowVector4d py_j,    // Y-positions of particle j
        double th_j,                // Thickness of particle j
        Eigen::RowVector3d do_j,    // Center of particle j
        Eigen::RowVector4d px_k,    // X-positions of particle k
        Eigen::RowVector4d py_k,    // Y-positions of particle k
        double th_k,                // Thickness of particle  k
        Eigen::RowVector3d do_k     // Center of particle k
        );

// Function to convert lenght width and distance to px py do_
void write_Points(
        double *w,
        double *l,
        double dx,
        Eigen::RowVector3d *d_or,
        Eigen::RowVector4d cortes_y,
        int j,
        Eigen::RowVector4d *px,
        Eigen::RowVector4d *py
        );

typedef Eigen::Matrix<double, 1, 7> RowVector7d; // Defines type to use in 's' Experiment of H_app

// Function to compute the Applied Field h_app based on 's'
Eigen::RowVector3d compute_Happ(
        int step,
        int n_of_particles
        );
#endif /* FUNCTIONS_SLLGS_H_ */
