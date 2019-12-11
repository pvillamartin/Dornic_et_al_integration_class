//
// Created by paula on 19/05/13.
//

#include <cstdio>
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include "Dornic_et_al_integration_method.h"

using namespace std;

int main(int argc, char *argv[]){
    ///////////////////////// DEFAULT DYNAMICS PARAMETERS ///////////////////////
    int number_cells=1024;                                      //Number of lattice cells
    double diffusion_coefficient=1;                             //Diffusion coefficient
    double linear_coefficient=0.01;                             //First order term coefficient
    double quadratic_coefficient=2;                             //Second order term coefficient
    double noise_coefficient=sqrt(1./double(number_cells));     //Demographic noise coefficient

    ///////////////////////// DEFAULT INTEGRATION PARAMETERS ///////////////////////
    double tmax=1000000;                                        //Maximum simulation time
    double dx=1;                                                //Spatial integration increment
    double dt=0.1;                                              //Temporal integration increment

    ///////////////////////// PARAMETERS FROM TERMINAL ///////////////////////
    // In order to alter the value of the parameters from terminal execute the
    // program in the following way:
    // ./main -D 1 -a 0.01 -b 2 -gamma 0.0001 -tmax 1000000 -dx 1 -dt 0.1
    for (int i = 1; i < argc; i++) {
        if (string(argv[i]) == "-D") {
            diffusion_coefficient = atof(argv[i + 1]);
        } else if (string(argv[i]) == "-a") {
            linear_coefficient = atof(argv[i + 1]);
        } else if (string(argv[i]) == "-b") {
            quadratic_coefficient = atof(argv[i + 1]);
        } else if (string(argv[i]) == "-gamma") {
            noise_coefficient = atof(argv[i + 1]);
        } else if (string(argv[i]) == "-tmax") {
            tmax = atof(argv[i + 1]);
        } else if (string(argv[i]) == "-dx") {
            dx = atof(argv[i + 1]);
        } else if (string(argv[i]) == "-dt") {
            dt = atof(argv[i + 1]);
        }
    }

    ///////////////////////// VARIABLES ///////////////////////
    vector <double> f_parameters;
    double t,mean_f;
    ofstream ft_file;
    char ft_name [100];
    sprintf(ft_name,"Dornic_integration_ft_D%g_a%g_b%g_tmax%g_dx%g_dt%g.txt",diffusion_coefficient,linear_coefficient,quadratic_coefficient,tmax,dx,dt);
    ft_file.open(ft_name, ios::out | ios::trunc);

    ///////////////////////// INITIALIZATION ///////////////////////
    //Random number generator
    int rand_seed=1;
    const gsl_rng_type * gsl_rng_T;
    gsl_rng * r;
    gsl_rng_env_setup();
    gsl_rng_default_seed=rand_seed;
    gsl_rng_T=gsl_rng_default;
    r=gsl_rng_alloc(gsl_rng_T);

    //Initialization of parameters
    f_parameters.push_back(diffusion_coefficient);
    f_parameters.push_back(linear_coefficient);
    f_parameters.push_back(noise_coefficient);
    f_parameters.push_back(quadratic_coefficient);

    //Declaration of Dornic class with the given parameters
    Dornic dornic(dt,dx,number_cells,&f_parameters);

    //Initialization of lattice
    for(int i=0;i<number_cells;i++) dornic.cell_density.push_back(1);
    dornic.set_1D_lattice();

    ///////////////////////// INTEGRATION ///////////////////////
    t=0;
    while(t<tmax){
        dornic.integrate(r);
        t+=dt;

        mean_f=0;
        for(int i=0;i<dornic.cell_density.size();i++) mean_f+=dornic.cell_density[i];
        ft_file<<t<<" "<<mean_f/number_cells<<endl;
    };

    return 0;
}

