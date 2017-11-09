//============================================================================
// Name        : s-LLGS.cpp
// Author      : João Nizer
// Version     : 0.0
// Copyright   : Universidade Federal de Minas Gerais - NanoComp
// Description : Solving s-LLGS in C++, Ansi-style
// References  : 2016 - Ament et al. - Solving the stochastic Landau-Lifshitz-
//			   : Gilbert-Slonczewski equation for monodomain nanomagnets.
//============================================================================

#include <iostream>
#include <stdio.h>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <ctime>
#include <chrono>
#include "matplotlibcpp.h"
#include <fstream>
#include <random>

#include "functions_sLLGS.h"
using namespace std;
namespace plt = matplotlibcpp;

//=========================== Constants ================================
double gama = 2.211E5; 			// Giromagnetic Ratio in 	[m/A.s]*
double q     = 1.60217662E-19; 	// Electron charge 			[C]
double mu0   = 4*atan(1)*4*1E-7;// [H/m] or 				[T.m/A]^
double kb    = 1.38064852E-23;	// Boltzman"s Constant 		[m2.kg/s2.K]
double hbar  = 1.0545718E-34;	// Norm. Planks Const. 		[J.s/rad]
double pi    = 4*atan(1);       // Pi 3.14
double t2am  = 1/mu0;			// Convert Tesla [T] to A/m
double j2ev  = 1/q;				// Convert Joule to eV
//* <cmath> already have a fun. called gamma with 2 M"s
//^ Pi value was calculated using atan(1)*4
//======================================================================

int main()
{

    int step=4000;                           // Time Steps
    int n_of_particles = 1 ; 				// Quantidade de particulas
    float Ms = 8E5;							// Saturation Manetization [A/m]
    float T = 0;							// Tmperature in Kelvin [K]
    double simulation_time = 40E-9; 		// Simulation Time in [s]
    double time_step = simulation_time/step;// Time step in [s]
    double dt = time_step*gama*Ms;			// Normalized Time_Step
    float alpha = 1; 						// Gilbert damping constant
    float alpha_l = 1/(1+pow(alpha,2));		// Scaling factor for LLGS
    double Ku = 500.0;                      // Anisotropy Constant
    clock_t begin,end;                      // Calculate Elapsed Time
    double elapsed_secs;                    // Elapsed Time
    cout << "Ms: " << Ms << "A/m" << endl;
    cout << "T: " << T << "k" << endl;
    cout << "Alpha" << alpha << endl;
    Eigen::RowVector3d n;					// Easy axis of particles
    Eigen::IOFormat formato(Eigen::FullPrecision); // Format Spec when printing Eigen vectors variable.format(formato)
    n << 0,1,0;


    Eigen::RowVector3d m[step][n_of_particles]; // Particles magnetization
    Eigen::RowVector3d h_app[step][n_of_particles]; // External applied field
    Eigen::RowVector3d i_s[step][n_of_particles]; // Spin current
    Eigen::Matrix3d     Nd      [n_of_particles];                   // Demagnetizing Tensor
    Eigen::Matrix3d     Nc      [n_of_particles][n_of_particles];   // Coupling Tensor
    Eigen::RowVector3d  h_c     [step][n_of_particles];             // Coupling field
    Eigen::RowVector3d  h_th    [step][n_of_particles];             // Thermal field
    Eigen::RowVector3d  h_eff   [step][n_of_particles];             // Effective field


    cout << "Cleaning Nd m h_eff h_app i_s h_c h_th" << endl;
    begin = clock();
    for (int i=0;i<step;i++){  					// Iterate over steps
        for (int j=0;j<n_of_particles;j++){ 	// Iterate over the particles in the step i
            for (int k=0;k<n_of_particles;k++){
                Nc[j][k]= Eigen::Matrix3d::Zero();
            }
            Nd[j]		= Eigen::Matrix3d::Zero();
            m[i][j] 	= Eigen::RowVector3d::Zero();//
            h_eff[i][j] = Eigen::RowVector3d::Zero();//
            h_app[i][j] = Eigen::RowVector3d::Zero();//
            i_s[i][j] 	= Eigen::RowVector3d::Zero();//
            h_c[i][j] 	= Eigen::RowVector3d::Zero();//
            h_th[i][j] 	= Eigen::RowVector3d::Zero();//
        }
    }

    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout << elapsed_secs << "s" << endl;

    cout << "Setting Initial Magnetization" << endl;
    begin = clock();
    for (int i=0;i<n_of_particles;i++){
        m[0][i] << 0     ,  pow(-1,i)*1.0000 ,         0; // Set initial magnetization for all particles;
    }

    cout << "m0.0: " << m[0][0] << " --- m0.1: " << m[0][1] << endl;

    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout << elapsed_secs << "s" << endl;

    //================ Specify the Applied Field

    int size_of_s = 4; // Number of specificaton lines used in 's'
    int sim_step_index;  // This variable acumulates the steps already used to describe the applied field
    // its maximum value is going to be equals to N -> total number of simulation steps

    double slope;


    RowVector7d s[n_of_particles][size_of_s];

    cout << "Setting H_app" << endl;
    begin = clock();
    //=========================================//
    //======== Experiment Description =========//
    //=========================================//
    //
    //  s[particle][line] << xi, yi, zi, xf, yf, zf, N
    // #i - initial value of h_app in Tesla on direction #
    // #f - final value of h_app in Tesla on direction #
    // N - number os sim steps the variation occurs
    //

    for (int i = 0; i <n_of_particles;i++){
        s[i][0] << 0, 0, 0, 0.25, 0, 0, step/4;
        s[i][1] << 0.25, 0, 0, 0.25, 0, 0, step/4;
        s[i][2] << 0.25, 0, 0, 0, 0, 0, step/4;
        s[i][3] << 0, 0, 0, 0, 0, 0, step/4;
    }



    //s[1][0] << 0, 0, 0, -1, -1, -1, step/4;
    //s[1][1] << -1, -1, -1, -1, -1, -1, step/4;
    //s[1][2] << -1, -1, -1, 0, 0, 0, step/4;
    //s[1][3] << 0, 0, 0, 0, 0, 0, step/4;

    //===== End of Experiment Description =====//



    for (int j=0; j<n_of_particles;j++){ // iterate on number of particles

        sim_step_index = 0; // Start index at 0

        for (int l=0; l< size_of_s; l++){ // iterate on lines of 's'

            for (int d = Coordinate::x ; d <= Coordinate::z; d++){ // iterate over the dimension (x,y,z)

                slope = (s[j][l](d+3)-s[j][l](d))/(s[j][l](6)-1); // Calculates the slope (x_f -x_i)/n;

                        for (int n=0; n < s[j][l](6); n++){ // iterate from 0 to N

                    h_app[sim_step_index+n][j](d) = (s[j][l](d) +n*slope)*t2am/Ms; // Caculates h_app = h_app_initial + n*slope
                    // and converts from T to A/m
                } //end iterate on N
            } // end iterate on dimensions
            sim_step_index += s[j][l](6); // Increment index base on the number the line 'l' lasts
            //cout << "inde: " << sim_step_index << endl;
        } // end iterate on lines of 's'
    } //end iterate on particles

    //for (int i=0; i<step;i++)
    //cout << i << ": " << h_app[i][0].format(formato) << endl;

    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout << elapsed_secs << "s" << endl;
    //================ Specify Particle Width, Length, Thickness and Distances in X-direction

    double w0 = 50;
    double l0 = 100;
    double th0 = 10;
    double dx0 = 10;

    double w [n_of_particles];
    double l [n_of_particles];
    double th[n_of_particles];
    double dx[n_of_particles];

    Eigen::RowVector4d px[n_of_particles];
    Eigen::RowVector4d py[n_of_particles];
    Eigen::RowVector4d cortes_y[n_of_particles];
    Eigen::RowVector3d do_[n_of_particles];



    //===============Set w, l, th, dx with defaults w0 l0 th0 and dx0

    cout << "Setting Particle Dimensions" << endl;
    begin = clock();

    for (int j=0; j<n_of_particles; j++){
        w[j] = w0;
        l[j] = l0;
        th[j] = th0;
        dx[j] = dx0;
        cortes_y[j] << 0,0,0,0;
    }


    for (int j=0; j<n_of_particles; j++){
        write_Points(w, l, dx[j], do_, cortes_y[j], j, px, py);
        //cout << do_[j] << endl;
        //cout << endl;
    }

    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout << elapsed_secs << "s" << endl;

    //================Calculate Demag Tensors

    cout << "Calculating Demag Tensors" << endl;
    begin = clock();

    for (int j=0;j<n_of_particles;j++){
        cout << j << " ";
        Nd[j]=compute_Demag( px[j], py[j], th[j]);
        for (int index=0; index<9; index++){
            Nd[j](index) = Nd[j](index)/(4*pi*w[j]*l[j]*th[j]);
        }
        Nd[j](3)=Nd[j](1);
        Nd[j](6)=Nd[j](2);
        Nd[j](7)=Nd[j](5);
        //cout << Nd[j].format(formato) << endl;
    }

    //    Nd[0]<<1.697951e-01,
    //            4.103516e-05,
    //            -2.310026e-02,
    //            4.103516e-05,
    //            8.306841e-02,
    //            1.265466e-02,
    //            -2.310026e-02,
    //            1.265466e-02,
    //            7.426134e-01;
    //    Nd[1]=Nd[0];

    cout << endl;
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout << elapsed_secs << "s" << endl;

    //================Compute Dipolar Tensors

    cout << "Calculating Dipolar Tensors" << endl;
    begin = clock();

    for (int j=0; j< n_of_particles; j++){
        // cout << j << " ";
        for (int k=j; k< n_of_particles; k++){
            if (j != k) {
                Nc[j][k]=compute_Dipolar( px[j], py[j], th[j], do_[j], px[k], py[k], th[k], do_[k]);
                for (int index=0;index<9;index++){
                    Nc[j][k](index) = Nc[j][k](index)/(4*pi*w[j]*l[j]*th[j]);
                }

                Nc[k][j]=Nc[j][k];

            }
            else {
                Nc[j][k] << 0,0,0,0,0,0,0,0,0;
            }

            // if (j==5 && k==10){
            //     cout << Nc[j][k].format(formato) << endl;
            //     cout << "j = " << j << endl;
            //     cout << "px: " << px[j] << "\n py: " << py[j] << "\n th: " << th[j] << "\n dor: " << do_[j] << endl;
            //     cout << "k = " << k << endl;
            //     cout << "px: " <<  px[k] << "\n py: " << py[k] << "\n th: " << th[k] << "\n dor: " << do_[k] << endl;
            //return 0;
            //}
            //cout << j << k << ":" << Nc[j][k].format(formato) << endl;
        }
    }
    cout << endl;
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout << elapsed_secs << "s" << endl;
    // Compute the Thermal Noise                      //
    // Standard Deviation
    double V = w[0]*l[0]*th[0]*1E-27; // If particle volume changes -> std must be reescaled

    double std = sqrt(2*alpha*kb*T/(gama*mu0*Ms*V*time_step));
    // Normalize the STD Deviation
    std = std/Ms;

    //
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution(0.0,std);
    cout << "std: " << std << "\nseed: " << seed << "\nrand sample: " << distribution(generator) << endl;
    for (int i=0;i < step;i++){
        for (int j=0;j < n_of_particles;j++){
            h_th[i][j] <<
                          distribution(generator),  // x
                    distribution(generator),  // y
                    distribution(generator); //z
        }
    }
    // Calculates the Anisotrpy Field Magnitude
    double H_k = 2*Ku/(Ms*mu0);
    double h_k = H_k / Ms; // Normalize H_k by Ms
    cout << h_k << endl;

    // Core Algorithm for Solving the s-LLGS equation //

    cout << "------------------------------------" 			<< endl;
    cout << "            Range-Kutta             " 			<< endl;
    cout << "------------------------------------" 			<< endl;
    cout << "dt real: \t\t" 	<< time_step 		<< " s"	<< endl;
    cout << "dt linha: \t\t" 	<< dt 		 				<< endl;
    cout << "tempo_total:\t" 	<< simulation_time 	<< " s" << endl;
    cout << "N: \t\t"           << step 					<< endl;
    cout << "------------------------------------"			<< endl;

    begin = clock();

    for (int i=0;i<step;i++){  					// Iterate over steps
        for (int j=0;j<n_of_particles;j++){ 	// Iterate over the particles in the step i

            for (int k=0;k<n_of_particles;k++){ // Calculates Coupling for step i
                h_c[i][j] += m[i][k]*Nc[j][k];
            } // end iterate over coupling

            // Compute the effective Field for particle j at step i
            h_eff[i][j] = compute_Heff(m[i][j], h_app[i][j], h_k, Nd[j],  h_c[i][j],  h_th[i][j], n);
            // Compute the new value of m for the step i+1;
            m[i+1][j] = rk4(m[i][j], h_eff[i][j], i_s[i][j], dt, alpha, alpha_l);

        } // end iterate over particles
    } // end iterate over steps
    // End of s-LLGS calculation
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout << elapsed_secs << "s" << endl;

    std::vector<double> time(step);
    std::vector<double> mx(step), my(step), mz(step);
    std::vector<double> hax(step), hay(step), haz(step);



    //    plt::figure();
    //    for (int j=0; j<n_of_particles;j++){

    //        for (int i=0;i<step;i++){
    //            int plot_particle = j; // chose particle to plot
    //            time[i] = i*time_step/1E-9;
    //            mx[i] = m[i][plot_particle](Coordinate::x);
    //            my[i] = m[i][plot_particle](Coordinate::y);
    //            mz[i] = m[i][plot_particle](Coordinate::z);
    //            hax[i] = h_app[i][plot_particle](Coordinate::x);
    //            hay[i] = h_c[i][plot_particle](Coordinate::y);
    //            haz[i] = h_c[i][plot_particle](Coordinate::z);
    //            //cout << mx[i] << " " << m[i][1](Coordinate::x) << endl;
    //        }
    //        cout << "m0: " << m[0][0] << " --- m_plot0: " << mx[0] << my[0] << mz[0] << endl;
    //        plt::subplot(n_of_particles,1,j+1);
    //        //plt::named_plot("m_y 1",time, mx, "b-");
    //        //plt::named_plot("m_y 2",time, my, "r-");
    //        //plt::named_plot("m_y 3",time, mz, "y-");
    //        //plt::named_plot("happ_x",time, hax, "b:");

    //        plt::plot(time, mx, "b-");
    //        plt::plot(time, my, "r-");
    //        plt::plot(time, mz, "y-");
    //        plt::plot(time, hax, "b:");

    //        plt::xlim(0.0, simulation_time/1E-9);
    //        plt::ylim(-1.2, 1.2);
    //        plt::ylabel("m");
    //        plt::legend();

    //        if (j==n_of_particles-1)
    //            plt::xlabel("Time (ns)");

    //        // Save the image (file format is determined by the extension)
    //        //if (j==0)

    //        //else
    //        //plt::save("../plots/particle1.png");
    //    }
    //    plt::save("../plots/particle0.png");


    // Export Data to File
    ofstream outputfile;
    outputfile.open("../plots/data.txt");

    outputfile << step << endl;
    outputfile << n_of_particles << endl;
    outputfile << time_step << endl;


    for (int i=0; i<step;i++){

        for (int j=0;j<n_of_particles;j++){

            outputfile << m[i][j](Coordinate::x) << "\t"  << m[i][j](Coordinate::y) << "\t" << m[i][j](Coordinate::z) << "\t" ;
        }
        outputfile << endl;
    }
}
