//
// Created by paula on 19/05/13.
// Updated by Victor on August 2020
//

#ifndef DORNIC_H
#define DORNIC_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#ifndef CONSTANT_CELL_NUMBER
#define CONSTANT_CELL_NUMBER true
#endif

#ifndef APPROX_POISSON
#define APPROX_POISSON 0
#endif

#ifndef RNG
#define RNG mt19937
#endif

#ifndef USE_EULER
#define USE_EULER false
#endif

using namespace std;

class Dornic
{

private:

    //RUNGE-KUTTA VARIABLES
    vector <double> k1,k2,k3,k4;
    vector <double> aux_cell_old, aux_cell_new;
    double dt,dx,dtm,dts;

    //DORNIC VARIABLES
    double lambda, lambda_product;
    poisson_distribution<int> poisson;
    gamma_distribution<double> gamma;
    normal_distribution<double> normal;

    //CELLS DENSITY AND ADJACENCY NETWORK
    vector <double> cell_density;
    double avg_density;
    int ncells;
    vector< vector<int> > neighbors;

    //COEFFICIENTS
    double linear_coeff, noise_coeff;

public:

    //Class constructor
    Dornic(const double dt_in, const double dx_in, const double cells_number, const vector<double> &f_parameters)
    {
        int i; //Counter

        //Initialization of integration increments
        dt=dt_in;
        dx=dx_in;

        //Initialization of Runge-Kutta variables
        dtm=0.5*dt;
        dts=dt/6.0;

        //Initial number of cells, set no activity for the system
        ncells = cells_number;
        cell_density = vector<double>(ncells, 0.0); 
        aux_cell_new = vector<double>(ncells);
        aux_cell_old = vector<double>(ncells);

        #if CONSTANT_CELL_NUMBER
        k1 = vector<double>(ncells, 0.0);
        k2 = vector<double>(ncells, 0.0);
        k3 = vector<double>(ncells, 0.0);
        k4 = vector<double>(ncells, 0.0);
        #endif

        //Initialization of Dornic variables
        set_coefficients(f_parameters);
    }

    //Overload with no parameters
    Dornic(const double dt_in, const double dx_in, const double cells_number)
    {
        int i; //Counter

        //Initialization of integration increments
        dt=dt_in;
        dx=dx_in;

        //Initialization of Runge-Kutta variables
        dtm=0.5*dt;
        dts=dt/6.0;

        //Initial number of cells, set no activity for the system
        ncells = cells_number;
        cell_density = vector<double>(ncells, 0.0); 
        aux_cell_new = vector<double>(ncells);
        aux_cell_old = vector<double>(ncells);

        #if CONSTANT_CELL_NUMBER
        k1 = vector<double>(ncells, 0.0);
        k2 = vector<double>(ncells, 0.0);
        k3 = vector<double>(ncells, 0.0);
        k4 = vector<double>(ncells, 0.0);
        #endif
    }



    ///////////////////////// DORNIC INTEGRATION ///////////////////////
    void integrate(RNG &gen)
    {
        //RUNGE-KUTTA VARIABLES
        #if !(CONSTANT_CELL_NUMBER)
        ncells = cell_density.size();
        k1 = vector<double>(ncells, 0.0);
        k2 = vector<double>(ncells, 0.0);
        k3 = vector<double>(ncells, 0.0);
        k4 = vector<double>(ncells, 0.0);
        #endif


        #if (!USE_EULER)
        //Runge-Kutta integration of the non-linear term and difussion.
        //Update of cells is done in the same loop as last RK step for efficiency
        f1(aux_cell_old, k1);
        f2f3(aux_cell_old, aux_cell_new, k2, dtm);
        aux_cell_old.swap(aux_cell_new);            //Swap contents is O(1), better than old = fast
        f2f3(aux_cell_old, aux_cell_new, k3, dt);
        aux_cell_old.swap(aux_cell_new);
        f4_and_stochastic(aux_cell_old, k1, k2, k3, gen);
        #else
        euler(aux_cell_old, gen);
        cell_density.swap(aux_cell_old); 
        #endif
    }

    //Non-constant coefficient overload of function
    void integrate(RNG &gen, const vector <double> &f_parameters)
    {
        //RUNGE-KUTTA VARIABLES
        #if !(CONSTANT_CELL_NUMBER)
        ncells = cell_density.size();
        f1 = vector<double>(ncells, 0.0);
        f2 = vector<double>(ncells, 0.0);
        f3 = vector<double>(ncells, 0.0);
        f4 = vector<double>(ncells, 0.0);
        #endif

        //DORNIC VARIABLES
        set_coefficients(f_parameters);

        #if (!USE_EULER)
        //Runge-Kutta integration of the non-linear term and difussion.
        //Update of cells is done in the same loop as last RK step for efficiency
        f1(aux_cell_old, k1);
        f2f3(aux_cell_old, aux_cell_new, k2, dtm);
        aux_cell_old.swap(aux_cell_new);            //Swap contents is O(1), better than old = fast
        f2f3(aux_cell_old, aux_cell_new, k3, dt);
        aux_cell_old.swap(aux_cell_new);
        f4_and_stochastic(aux_cell_old, k1, k2, k3, gen);
        #else
        euler(aux_cell_old, gen);
        cell_density.swap(aux_cell_old); 
        #endif
    }

    ///////////////////////// INITIALIZATION FUNCTIONS ///////////////////////

    //Set all the parameters
    void set_coefficients(const vector <double> &f_parameters)
    {
        set_essential_coefficients(f_parameters);
        set_non_linear_coefficients(f_parameters);

        double lambda_const=2./(noise_coeff*noise_coeff);
        double lambda_exp=exp(-linear_coeff*dt);
        lambda=lambda_const*linear_coeff*lambda_exp/(1.0-lambda_exp);
        lambda_product=lambda/lambda_exp;
    }

    //Declaration, to be filled by user
    void set_non_linear_coefficients(const vector <double> &f_parameters);

    //Coefficients needed to integrate linear + noise parts
    void set_essential_coefficients(const vector <double> &f_parameters)
    {
        linear_coeff = f_parameters[0];
        noise_coeff = f_parameters[1];

        double lambda_const=2./(noise_coeff*noise_coeff);
        double lambda_exp=exp(-linear_coeff*dt);
        lambda=lambda_const*linear_coeff*lambda_exp/(1.0-lambda_exp);
        lambda_product=lambda/lambda_exp;
    }


    ///////////////////////// CORE INTEGRATION FUNCTIONS ///////////////////////

    //Declaration -user will fill this function in the CPP file
    double non_linear_rhs(const int inode, const vector<double> &field) const;

    //RK first function
    void f1(vector<double> &aux_cell, vector<double> &k1)
    {
        int i;

        for(i=0; i<ncells; i++)
        {
            k1[i] = non_linear_rhs(i, cell_density);
            aux_cell[i] = cell_density[i] + dtm*k1[i];
        }
    }

    //RK second and third functions
    void f2f3(const vector<double> &aux_old, vector<double> &aux_new, vector<double> &k_out, const double dt_in)
    {
        int i;

        for(i=0; i<ncells; i++)
        {
            k_out[i] = non_linear_rhs(i, aux_old);
            aux_new[i] = cell_density[i] + dt_in*k_out[i];
        }
    }

    //RK fourth function and stochastic step, all in the same loop
    void f4_and_stochastic(const vector<double> &aux_old, const vector<double> &k1, const vector<double> &k2, const vector<double> &k3, RNG &gen)
    {
        int i;

        double k4, cell_new;
        double mu;

        avg_density = 0.0;
        for(i=0; i<ncells; i++)
        {
            k4 = non_linear_rhs(i, aux_old);

            cell_density[i] += dts * (k1[i] + 2*(k2[i] + k3[i]) + k4);

            #if APPROX_POISSON==0
            poisson = poisson_distribution<int>(lambda_product * cell_density[i]);
            gamma = gamma_distribution<double>(poisson(gen), 1.0);
            #else
            mu = lambda_product * cell_density[i];
            if (mu > APPROX_POISSON)
            {
                normal = normal_distribution<double>(mu, sqrt(mu));
                gamma = gamma_distribution<double>(normal(gen), 1.0);
            }
            else
            {
                poisson = poisson_distribution<int>(lambda_product * cell_density[i]);
                gamma = gamma_distribution<double>(poisson(gen), 1.0);
            }
            #endif
            cell_density[i]= gamma(gen)/lambda;

            //Make averages
            avg_density += cell_density[i];
        }    
        avg_density /= 1.0 * ncells;    
    }

    //Euler integration
    void euler(vector<double> &aux, RNG &gen)
    {
        int i;
        double mu;
        double f;

        avg_density = 0.0;
        for(i=0; i<ncells; i++)
        {
            f = non_linear_rhs(i, cell_density);

            aux[i] = cell_density[i] + dt * f;

            #if APPROX_POISSON==0
            poisson = poisson_distribution<int>(lambda_product * aux[i]);
            gamma = gamma_distribution<double>(poisson(gen), 1.0);
            #else
            mu = lambda_product * aux[i];
            if (mu > APPROX_POISSON)
            {
                normal = normal_distribution<double>(mu, sqrt(mu));
                gamma = gamma_distribution<double>(normal(gen), 1.0);
            }
            else
            {
                poisson = poisson_distribution<int>(lambda_product * aux[i]);
                gamma = gamma_distribution<double>(poisson(gen), 1.0);
            }
            #endif

            aux[i]= gamma(gen)/lambda;

            //Make averages
            avg_density += aux[i];
        }    
        avg_density /= 1.0 * ncells;   
    }

    ///////////////////////// NETWORK ///////////////////////
    void set_1D_lattice(const bool periodic = false)
    {
        int i;

        //Initialize vector with size 2
        neighbors = vector<vector<int>>(ncells, vector<int>(2));

        //Bulk, set all neighbours
        for(i=1; i<ncells-1; i++)
        {
            neighbors[i][0] = i-1;
            neighbors[i][1] = i+1;
        }

        //At the end of the lattice, set the boundaries
        if (not periodic)
        {
            //These ones only have size one, re-define
            neighbors[0] = vector<int>(1, 1);
            neighbors[ncells-1] = vector<int>(1, ncells-2);
        }
        else
        {
            neighbors[0][0] = ncells-1;
            neighbors[0][1] = 1;
            neighbors[ncells-1][1] = ncells-2;
            neighbors[ncells-1][1] = 0;        
        }
    }

    void set_2D_lattice(const bool periodic = false)
    {
        double N=cell_density.size();

        int x,y,i; //Counters
        int up, down, right, left;
        int L = sqrt(cell_density.size());

        neighbors = vector<vector<int>>(ncells, vector<int>(4));

        //Set all neighbours
        if (not periodic)
        {
            for(y=1; y < L-1; y++)
            {
                for (x=1; x < L-1; x++)
                {
                    i = x+y*L;

                    neighbors[i][0] = x + (y+1)*L; //Up
                    neighbors[i][1] = x + (y-1)*L; //Down
                    neighbors[i][2] = (x+1) + y*L; //Left
                    neighbors[i][3] = (x-1) + y*L; //Right
                }
            }

            for (i=1; i < L-1; i++)
            {
                down = i; //Bottom row
                up = i + (L-1)*L; //Top row
                left = i*L; //Left column
                right = (L-1) + i*L; //Right column

                //Each one of these have three neighbours, redefine
                neighbors[down] = vector<int>(3);
                neighbors[up] = vector<int>(3);
                neighbors[left] = vector<int>(3);
                neighbors[right] = vector<int>(3);

                //Set the neighbours
                neighbors[down][0] = down + L;
                neighbors[down][1] = down - 1; 
                neighbors[down][2] = down + 1; 

                neighbors[up][0] = up - L;
                neighbors[up][1] = up - 1; 
                neighbors[up][2] = up + 1;

                neighbors[left][0] = left + L;
                neighbors[left][1] = left - L; 
                neighbors[left][2] = left + 1;

                neighbors[right][0] = right + L;
                neighbors[right][1] = right - L; 
                neighbors[right][2] = right - 1;
            }

            //Finally, make the corners
            neighbors[0][0] = 1; 
            neighbors[0][1] = L;

            neighbors[L-1][0] = L-2; 
            neighbors[L-1][1] = (L-1)+L;

            neighbors[(L-1)*L][0] = 1+(L-1)*L; 
            neighbors[(L-1)*L][1] = (L-2)*L;

            neighbors[(L-1)+(L-1)*L][0] = (L-2)+(L-1)*L; 
            neighbors[(L-1)+(L-1)*L][1] = (L-1)+(L-2)*L;
        }
        else
        {

            //A bit slower than the code above, but constructs the full system
            for(y=0; y < L; y++)
            {
                for (x=0; x < L; x++)
                {
                    i = x+y*L;

                    up = (y < L-1) ? x + (y+1)*L : x;
                    down = (y > 0) ? x + (y-1)*L : x + (L-1)*L;
                    right = (x < L-1) ? x+1 + y*L : y*L;
                    left = (x > 0) ? x-1 + y*L : L-1 + y*L;

                    neighbors[i][0] = up; //Up
                    neighbors[i][1] = down; //Down
                    neighbors[i][2] = right; //Left
                    neighbors[i][3] = left; //Right
                }
            }
        }

    }

    void custom_network(const vector<vector<int> > &network)
    {
        neighbors = network;
    }

    ///////////////////////// INITIAL CONDITIONS ///////////////////////

    //Random initial values for all cells between minv and maxv
    void random_intial_cond(RNG &gen, const double minv = 0.0, const double maxv = 1.0)
    {
        uniform_real_distribution<double> unif(minv, maxv);
        int i;

        avg_density = 0.0;
        for (i=0; i < cell_density.size(); i++)
        {
            cell_density[i] += unif(gen);
            avg_density += cell_density[i];
        }
        avg_density /= 1.0*ncells;
    }

    //All cells same value
    void homogeneous_initial_cond(const double density_value = 1.0)
    {
        cell_density = vector<double>(ncells, density_value);
        avg_density = density_value;
    }

    //All cells to 0 except the specified one
    void single_seed(const int inode, const double value = 1.0)
    {
        cell_density[inode] = value;
        avg_density = value / (1.0 * ncells);
    }    

    ///////////////////////// OPERATOR OVERLOADS ///////////////////////

    double &operator[](const int index)
    {
        return cell_density[index];
    }

    //The same as above, but this will be invoking when reading the value.
    const double &operator[](const int index) const 
    {
        return cell_density[index];
    }

    //Declare as friend function
    friend ostream& operator<<(ostream& os, const Dornic &dornic);

    ///////////////////////// USABILITY ///////////////////////

    //Make visible the density
    double density() 
    {
        return avg_density;
    }

    //Check the average mean of the Poisson distribution
    double avg_poisson_mean()
    {
        return lambda_product * avg_density;
    }

};

#endif