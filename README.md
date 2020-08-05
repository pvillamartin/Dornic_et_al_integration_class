# Dornic integration method for multipicative noise

***********************************************************************************************************************************
## Overview

`Dornic_et_al_integration_method.h` provides a general class to integrate stochastic differential equations.

It is an optimization and generalization of a basic code implemented by Ivan Dornic and Juan A. Bonachela based on the algorithm proposed in:

> “Integration of langevin equations with multiplicative noise and the viability of field theories for absorbing  phase transitions”; Ivan Dornic, Hugues Chaté, and Miguel A Munoz;  Physical review letters 94; 100601 (2005).

and improved following:

> “Simulation of spatial systems with demographic noise”; Haim Weissmann, Nadav M Shnerb, and David A Kessler; Physical Review E 98; 022131 (2018).

This project is licensed under the terms of the GNU General Public License license. If you use our code in an academic publication, please include us in the acknowledgements and cite the previous papers.

***********************************************************************************************************************************

## Installation and compiling

Download or clone this repository using `git clone`. Integration routines are inside the file`Dornic_et_al_integration_method.h`, that should included in your C++ file. Then compile the code as usual, 

```bash
g++ example_main.cpp -std=c++11 -other_flags -o example_main
```
The Dornic integration method provides also two compilation options: `CONSTANT_CELL_NUMBER` and `RNG`. Their default values are `true` and `mt19937` by default. The `CONSTANT_CELL_NUMBER` macro states if the algorithm parameters are constant during integration time or are functions of time itself. If you intend to change parameters during integration, set it to `false`, as

```bash
g++ example_main.cpp -std=c++11 -other_flags -DCONSTANT_CELL_NUMBER=false -o example_main
```

Finally, the `RNG` is a macro that defines the currently used random number generator. Any valid C++ generator can be set on compilation. 

***********************************************************************************************************************************

##  Use

In order to show the proper way of using it, we implemented the `example_main.cpp`. These codes integrate  the stochastic equation described in "Dornic_method.pdf" in a 1D lattice.  

1. To integrate a stochastic equation with different non-linear terms you may need to change the functions:

    ```c++
    ///////////////////////// NON-LINEAR TERM INTEGRATION FUNCTIONS ///////////////////////
    	void set_non_linear_coefficients(const vector <double> &f_parameters){
    		quadratic_coefficient=f_parameters[3];
    	}
    	double non_linear_term_integration(const int inode, const vector <double> &f_in, const double dt_in){
    		return -quadratic_coefficient*(cell_density[inode]+dt_in*(*f_in)[inode])*(cell_density[inode]+dt_in*(*f_in)[inode]);
    	}
    ```

    

2. To integrate a stochastic equation with a different adyacency network, create a new function similar to the following:

    ```c++
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
    ```

    

3. To integrate a differential equation with coefficients that change in time, just use the method `f_dornic.integrate(r,f_parameters);`



4. To integrate a differential equation with number of cells that change in time, make `#define CONSTANT_CELL_NUMBER` to `false` as indicated above.


***********************************************************************************************************************************
