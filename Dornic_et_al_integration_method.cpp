//
// Created by paula on 19/05/13.
// Updated by Victor on August 2020
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>

using namespace std;

#include "Dornic_et_al_integration_method.h"

//User-defined variables
double quad_coeff;
double D;

void Dornic::set_non_linear_coefficients(const vector <double> &f_parameters)
{
    quad_coeff = f_parameters[2];
    D = f_parameters[3] / (dx * dx);
}

double Dornic::non_linear_rhs(const int inode, const vector<double> &field) const
{
    //Non-linear terms
    const double quadratic = -quad_coeff * field[inode] * field[inode];

    //Diffusion
    double diff_sum;
    int num_neigh, index_neigh;

    num_neigh = neighbors[inode].size();
    diff_sum = 0.0;

    //DIFFUSION INTEGRATION
    for(int i=0; i<num_neigh; i++)
    {
        index_neigh = neighbors[inode][i];
        diff_sum += field[index_neigh];
    }

    diff_sum = D*(diff_sum - num_neigh*field[inode]);

    return diff_sum + quadratic;
}

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
