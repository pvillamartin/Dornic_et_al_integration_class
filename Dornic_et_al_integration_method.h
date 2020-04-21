//
// Created by paula on 19/05/13.
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
#include<random>

#define CONSTANT_COEFFICIENTS true
#define CONSTANT_CELLS_NUMBER true
#define RNG mt19937

using namespace std;

class Dornic{

public:
    //CELLS DENSITY AND ADJACENCY NETWORK
    vector <double> cell_density;
    struct nodes_struct{
        vector<int> neighbors;
    };
    vector <nodes_struct> nodes_neighbors_vector;
    //COEFFICIENTS
    double diffusion_coefficient,linear_coefficient,noise_coefficient,quadratic_coefficient;
    //RUNGE-KUTTA VARIABLES
    vector <double> f1,f2,f3,f4;
    double dt,dx,dtm,dts;
    //DORNIC VARIABLES
    double lambda,lambda_product;
    //Class constructor
    Dornic(double dt_in, double dx_in,double cells_number,vector <double> *f_parameters=NULL){

        //Initialization of integration increments
        dt=dt_in;
        dx=dx_in;

        //Initialization of the adjacency network
        for(int i=0;i<cells_number;i++){
            nodes_struct empty_neighbors_vector;
            nodes_neighbors_vector.push_back(empty_neighbors_vector);
        }

        //Initialization of Runge-Kutta variables
        dtm=0.5*dt;
        dts=dt/6.0;
        #if CONSTANT_CELLS_NUMBER

        f1 = vector<double>(cells_number, 0.0);
        f2 = vector<double>(cells_number, 0.0);
        f3 = vector<double>(cells_number, 0.0);
        f4 = vector<double>(cells_number, 0.0);
        #endif

        //Initialization of Dornic variables
        #if CONSTANT_COEFFICIENTS
        set_coefficients(f_parameters);
        double lambda_const=2./(noise_coefficient*noise_coefficient);
        double lambda_exp=exp(linear_coefficient*dt);
        lambda=lambda_const*linear_coefficient/(lambda_exp-1.);
        lambda_product=lambda*lambda_exp;
        #endif
    }

    ///////////////////////// DORNIC INTEGRATION ///////////////////////
    void integrate(RNG &gen, vector <double> *f_parameters=NULL){
        //RUNGE-KUTTA VARIABLES
        #if !(CONSTANT_CELLS_NUMBER)
        f1.clear();
        f2.clear();
        f2.clear();
        f2.clear();

        f1 = vector<double>(cell_density->size(), 0.0);
        f2 = vector<double>(cell_density->size(), 0.0);
        f3 = vector<double>(cell_density->size(), 0.0);
        f4 = vector<double>(cell_density->size(), 0.0);
        #endif


        //DORNIC VARIABLES
        #if !(CONSTANT_COEFFICIENTS)
        set_coefficients(f_parameters);
        double lambda_const=2./(noise_coefficient*noise_coefficient);
        double lambda_exp=exp(linear_coefficient*dt);
        lambda=lambda_const*linear_coefficient/(lambda_exp-1.);
        lambda_product=lambda*lambda_exp;
        #endif

        ///Runge-Kutta integration of the non-linear term and difussion
        RungeKutta_integrate(&f1, &f1, 0);
        RungeKutta_integrate(&f1, &f2, dtm);
        RungeKutta_integrate(&f2, &f3, dtm);
        RungeKutta_integrate(&f3, &f4, dt);

        //Define distributions to be used later
        poisson_distribution<int> poisson;
        gamma_distribution<double> gamma;

        ///Dornic integration of the noise and linear term
        for(int i=0;i<cell_density.size();i++){
            cell_density[i]+=dts*(f1[i]+2.0*f2[i]+2.0*f3[i]+f4[i]);

            poisson = poisson_distribution<int>(lambda_product * cell_density[i]);
            gamma = gamma_distribution<double>(poisson(gen), 1.0);

            cell_density[i]= gamma(gen)/lambda;
        }
    }

    ///////////////////////// BASIC INTEGRATION FUNCTIONS ///////////////////////
    void set_coefficients(vector <double> *f_parameters){
        set_essential_coefficients(f_parameters);
        set_non_linear_coefficients(f_parameters);
    }
    void set_essential_coefficients(vector <double> *f_parameters){
        diffusion_coefficient=(*f_parameters)[0];
        linear_coefficient=(*f_parameters)[1];
        noise_coefficient=(*f_parameters)[2];
    }
    double diffusion_term_integration(int inode, vector <double> *f_in, double dt_in){
        //PARAMETERS
        double D=diffusion_coefficient/(dx*dx);

        //VARIABLES
        double diffusion_integrated_term,diffusion_sum_f=0,diffusion_sum_fin=0;
        int neighbors=0,ineighbor;

        //DIFFUSION INTEGRATION
        for(int i=0;i<nodes_neighbors_vector[inode].neighbors.size();i++){
            ineighbor=nodes_neighbors_vector[inode].neighbors[i];
            diffusion_sum_f+=cell_density[ineighbor];
            diffusion_sum_fin+=(*f_in)[ineighbor];
            neighbors++;
        }
        diffusion_integrated_term=D*(diffusion_sum_f-neighbors*cell_density[inode]+dt_in*(diffusion_sum_fin-neighbors*(*f_in)[inode]));
        return diffusion_integrated_term;
    }
    void RungeKutta_integrate(vector<double> *f_in,vector<double> *f_out, double dt_in){
        for(int i=0;i<cell_density.size();i++){
            (*f_out)[i]=diffusion_term_integration(i,f_in,dt_in)+non_linear_term_integration(i,f_in,dt_in);
        }
    }

    ///////////////////////// NON-LINEAR TERM INTEGRATION FUNCTIONS ///////////////////////
    void set_non_linear_coefficients(vector <double> *f_parameters){
        quadratic_coefficient=(*f_parameters)[3];
    }
    double non_linear_term_integration(int inode, vector <double> *f_in, double dt_in){
        return -quadratic_coefficient*(cell_density[inode]+dt_in*(*f_in)[inode])*(cell_density[inode]+dt_in*(*f_in)[inode]);
    }

    ///////////////////////// NETWORK ///////////////////////
    void set_1D_lattice(){
        double N=cell_density.size();
        for(int inode=0;inode<N;inode++) nodes_neighbors_vector[inode].neighbors.clear();
        for(int inode=1;inode<N-1;inode++){
            nodes_neighbors_vector[inode].neighbors.push_back(inode-1);
            nodes_neighbors_vector[inode].neighbors.push_back(inode+1);
        }

        nodes_neighbors_vector[0].neighbors.push_back(1);
        nodes_neighbors_vector[N-1].neighbors.push_back(N-2);
    }
};






