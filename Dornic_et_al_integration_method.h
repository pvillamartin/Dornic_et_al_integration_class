//
// Created by paula on 19/05/13.
// Updated by Victor on August 2020
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
#include <random>

#ifndef CONSTANT_CELL_NUMBER
#define CONSTANT_CELL_NUMBER true
#endif

#ifndef RNG
#define RNG mt19937
#endif

using namespace std;

class Dornic{

public:
    //CELLS DENSITY AND ADJACENCY NETWORK
    vector <double> cell_density;
    double avg_density;
    int ncells;
    vector< vector<int> > neighbors;

    //COEFFICIENTS
    double D, linear_coeff, noise_coeff, quadratic_coeff;

    //RUNGE-KUTTA VARIABLES
    vector <double> f1,f2,f3,f4;
    double dt,dx,dtm,dts;

    //DORNIC VARIABLES
    double lambda,lambda_product;

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

        #if CONSTANT_CELL_NUMBER
        f1 = vector<double>(ncells, 0.0);
        f2 = vector<double>(ncells, 0.0);
        f3 = vector<double>(ncells, 0.0);
        f4 = vector<double>(ncells, 0.0);
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

        #if CONSTANT_CELL_NUMBER
        f1 = vector<double>(ncells, 0.0);
        f2 = vector<double>(ncells, 0.0);
        f3 = vector<double>(ncells, 0.0);
        f4 = vector<double>(ncells, 0.0);
        #endif
    }



    ///////////////////////// DORNIC INTEGRATION ///////////////////////
    void integrate(RNG &gen)
    {
        //RUNGE-KUTTA VARIABLES
        #if !(CONSTANT_CELL_NUMBER)
        ncells = cell_density.size();
        f1 = vector<double>(ncells, 0.0);
        f2 = vector<double>(ncells, 0.0);
        f3 = vector<double>(ncells, 0.0);
        f4 = vector<double>(ncells, 0.0);
        #endif

        ///Runge-Kutta integration of the non-linear term and difussion
        RungeKutta_integrate(f1, f1, 0);
        RungeKutta_integrate(f1, f2, dtm);
        RungeKutta_integrate(f2, f3, dtm);
        RungeKutta_integrate(f3, f4, dt);

        //Define distributions to be used later
        poisson_distribution<int> poisson;
        gamma_distribution<double> gamma;

        ///Dornic integration of the noise and linear term
        avg_density = 0.0;

        for(int i=0; i<ncells; i++)
        {
            cell_density[i]+=dts*(f1[i]+2.0*f2[i]+2.0*f3[i]+f4[i]);

            poisson = poisson_distribution<int>(lambda_product * cell_density[i]);
            gamma = gamma_distribution<double>(poisson(gen), 1.0);

            cell_density[i]= gamma(gen)/lambda;

            //Make averages
            avg_density += cell_density[i];
        }

        //Finish averages
        avg_density /= 1.0*ncells;
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

        ///Runge-Kutta integration of the non-linear term and difussion
        RungeKutta_integrate(f1, f1, 0);
        RungeKutta_integrate(f1, f2, dtm);
        RungeKutta_integrate(f2, f3, dtm);
        RungeKutta_integrate(f3, f4, dt);

        //Define distributions to be used later
        poisson_distribution<int> poisson;
        gamma_distribution<double> gamma;

        ///Dornic integration of the noise and linear term
        avg_density = 0.0;

        for(int i=0; i<ncells; i++)
        {
            cell_density[i]+=dts*(f1[i]+2.0*f2[i]+2.0*f3[i]+f4[i]);

            poisson = poisson_distribution<int>(lambda_product * cell_density[i]);
            gamma = gamma_distribution<double>(poisson(gen), 1.0);

            cell_density[i]= gamma(gen)/lambda;


            //Make averages
            avg_density += cell_density[i];
        }
        //Finish averages
        avg_density /= 1.0*ncells;
    }


    ///////////////////////// BASIC INTEGRATION FUNCTIONS ///////////////////////
    void set_coefficients(const vector <double> &f_parameters)
    {
        set_essential_coefficients(f_parameters);
        set_non_linear_coefficients(f_parameters);

        double lambda_const=2./(noise_coeff*noise_coeff);
        double lambda_exp=exp(-linear_coeff*dt);
        lambda=lambda_const*linear_coeff*lambda_exp/(1.0-lambda_exp);
        lambda_product=lambda/lambda_exp;

        cout << lambda_const << " " << lambda_exp << " " << lambda << " " << lambda_product << endl;
    }


    void set_essential_coefficients(const vector <double> &f_parameters)
    {
        D = f_parameters[0] / (dx * dx);
        linear_coeff = f_parameters[1];
        noise_coeff = f_parameters[2];
    }

    double diffusion_term_integration(const int inode, const vector <double> &f_in, const double dt_in)
    {
        //PARAMETERS
        //double D=diffusion_coefficient/(dx*dx);

        //VARIABLES
        double diff_integrated, dif_sum_f, dif_sum_fin;
        int num_neigh, index_neigh;

        num_neigh = neighbors[inode].size();
        dif_sum_f = dif_sum_fin = 0.0;

        //DIFFUSION INTEGRATION
        for(int i=0; i<num_neigh; i++)
        {
            index_neigh= neighbors[inode][i];
            dif_sum_f += cell_density[index_neigh];
            dif_sum_fin += f_in[index_neigh];
        }
        
        diff_integrated = D*(dif_sum_f - num_neigh*cell_density[inode] + dt_in * (dif_sum_fin - num_neigh * f_in[inode]));
        return diff_integrated;
    }

    void RungeKutta_integrate(const vector<double> &f_in, vector<double> &f_out, const double dt_in)
    {
        int i;

        for(i=0; i<ncells; i++)
        {
            f_out[i] = diffusion_term_integration(i,f_in,dt_in) + non_linear_term_integration(i,f_in,dt_in);
        }
    }

    ///////////////////////// NON-LINEAR TERM INTEGRATION FUNCTIONS ///////////////////////
    void set_non_linear_coefficients(const vector <double> &f_parameters)
    {
        quadratic_coeff=f_parameters[3];
    }


    double non_linear_term_integration(const int inode, const vector <double> &f_in, const double dt_in)
    {
        double auxiliar = (cell_density[inode] + dt_in*f_in[inode]);
        return -quadratic_coeff * auxiliar * auxiliar;
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

    void homogeneous_initial_cond(const double density_value = 1.0)
    {
        cell_density = vector<double>(ncells, density_value);
        avg_density = density_value;
    }

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
};

//Ostream operator allows to write all the object information
ostream& operator<<(ostream& os, const Dornic &dornic)
{
    int i;
    for (i=0; i < dornic.ncells; i++)
    {
        os << dornic[i] << " ";
    }
    os << endl;
    return os;
}