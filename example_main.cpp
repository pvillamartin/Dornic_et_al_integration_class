//
// Created by paula on 19/05/13.
//

#include <cstdio>
#include <iostream>
#include <vector>
#include <random>
#include "Dornic_et_al_integration_method.h"

using namespace std;

int main(int argc, char *argv[])
{

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
    for (int i = 1; i < argc; i++) 
    {
        if (string(argv[i]) == "-D") diffusion_coefficient = atof(argv[i + 1]);
        else if (string(argv[i]) == "-a") linear_coefficient = atof(argv[i + 1]); 
        else if (string(argv[i]) == "-b") quadratic_coefficient = atof(argv[i + 1]);
        else if (string(argv[i]) == "-gamma") noise_coefficient = atof(argv[i + 1]);    
        else if (string(argv[i]) == "-tmax") tmax = atof(argv[i + 1]); 
        else if (string(argv[i]) == "-dx") dx = atof(argv[i + 1]);
        else if (string(argv[i]) == "-dt") dt = atof(argv[i + 1]);
    }

    ///////////////////////// VARIABLES ///////////////////////

    //Random number generator
    int rand_seed=45345456;
    RNG gen(rand_seed); //MT19937 is a RNG with long period, good quality

    //Init variables
    vector <double> f_parameters(4);
    double t, mean_f;
    ofstream ft_file;

    ///////////////////////// INITIALIZATION ///////////////////////

    //Set parameters
    f_parameters = {diffusion_coefficient, linear_coefficient, noise_coefficient, quadratic_coefficient};

    //Declaration of Dornic class with the given parameters
    Dornic dornic(dt,dx,number_cells,f_parameters);

    //Initialization of lattice and initial conditions
    dornic.set_1D_lattice();
    dornic.random_intial_cond(gen);

    ///////////////////////// INTEGRATION ///////////////////////

    ft_file.open("prueba");
    for (t=0; t<=tmax; t+=dt)
    {
        dornic.integrate(gen);
        t+=dt;
        ft_file << t << " " << dornic.avg_density <<endl;
    }
    ft_file.close();

    return EXIT_SUCCESS;
}

