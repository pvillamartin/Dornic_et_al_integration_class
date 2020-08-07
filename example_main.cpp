//
// Created by paula on 19/05/13.
//

#include <cstdio>
#include <iostream>
#include <vector>
#include <random>
#include <cstring>
#include "Dornic_et_al_integration_method.h"

using namespace std;

#ifndef DIAGRAM
#define DIAGRAM false
#endif


int main(int argc, char *argv[])
{

    ///////////////////////// DEFAULT DYNAMICS PARAMETERS ///////////////////////
    int number_cells=4096;                                      //Number of lattice cells
    double diffusion_coefficient=0.1;                           //Diffusion coefficient
    double linear_coefficient=0.01;                             //First order term coefficient
    double quadratic_coefficient=2;                             //Second order term coefficient
    double noise_coefficient=1.0;                               //Demographic noise coefficient

    ///////////////////////// DEFAULT INTEGRATION PARAMETERS ///////////////////////
    double tmax=100;                                        //Maximum simulation time
    double dx=0.5;                                          //Spatial integration increment
    double dt=0.01;                                         //Temporal integration increment

    ///////////////////////// OTHER PARAMETERS   ///////////////////////////////////
    string filepath = "output";                             //Where to write results

    ///////////////////////// PARAMETERS FROM TERMINAL ///////////////////////
    // In order to alter the value of the parameters from terminal execute the
    // program in the following way:
    // ./main -D 1 -a 0.01 -b 2 -gamma 0.0001 -tmax 1000000 -dx 1 -dt 0.1 -path myfile.txt
    for (int i = 1; i < argc; i++) 
    {
        if (string(argv[i]) == "-D") diffusion_coefficient = atof(argv[i + 1]);
        else if (string(argv[i]) == "-a") linear_coefficient = atof(argv[i + 1]); 
        else if (string(argv[i]) == "-b") quadratic_coefficient = atof(argv[i + 1]);
        else if (string(argv[i]) == "-gamma") noise_coefficient = atof(argv[i + 1]);    
        else if (string(argv[i]) == "-tmax") tmax = atof(argv[i + 1]); 
        else if (string(argv[i]) == "-dx") dx = atof(argv[i + 1]);
        else if (string(argv[i]) == "-dt") dt = atof(argv[i + 1]);
        else if (string(argv[i]) == "-path") filepath = string(argv[i + 1]);

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
    f_parameters = {linear_coefficient, noise_coefficient, quadratic_coefficient, diffusion_coefficient};

    //Declaration of Dornic class with the given parameters
    Dornic dornic(dt,dx,number_cells);

    //Initialization of lattice and initial conditions
    dornic.set_2D_lattice();

    ///////////////////////// INTEGRATION ///////////////////////

    #if !(DIAGRAM)

    dornic.set_coefficients(f_parameters);
    dornic.random_intial_cond(gen);

    tmax = 500.0;
    
    ft_file.open(filepath);
    for (t=0; t<=tmax; t+=dt)
    {
        dornic.integrate(gen);
        ft_file << t << " " << dornic.density() <<endl;
    }
    ft_file.close();

    #else

    double lambda;          //Parameter
    double aver, aver2;     //To make averages
    double t_heatup = 10.0; //Time to thermalize
    double ro_total;        //Total activity in the system
    double thrs = 1e-2;     //For speedup, terminate if activity too low 

    int nits;
    ft_file.open(filepath);
    //lambda = 1.1;
    for (lambda = 1.1; lambda <= 2.0; lambda += 0.05)
    {
        cout << lambda << endl;

        //Contact-Process coefficients
        linear_coefficient = lambda - 1.0;
        quadratic_coefficient = lambda;

        //Re-set coefficients and re-start integration
        f_parameters[0] = linear_coefficient;
        f_parameters[2] = quadratic_coefficient;

        dornic.set_coefficients(f_parameters);
        dornic.random_intial_cond(gen);

        //Integrate and make averages
        aver = aver2 = 0.0;
        t = 0;
        ro_total = 1.0;
        nits = 0;
        while (t < tmax && ro_total > thrs)
        {
            dornic.integrate(gen);
            if (t > t_heatup)
            {
                aver += dornic.density();
                aver2 += dornic.density() * dornic.density();
                nits++;
            }
            ro_total = dornic.density() * number_cells;

            t += dt;
        }

        if (ro_total > thrs)
        {
            aver /= nits * 1.0;
            aver2 /= nits * 1.0;
        }
        else
        {
            aver = aver2 = 0.0;
        }
        

        ft_file << lambda << " " << aver << " " << aver2 - aver*aver << endl;
    }
    ft_file.close();
    #endif

    return EXIT_SUCCESS;
} 

