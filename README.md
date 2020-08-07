# Dornic integration method for multipicative noise

***********************************************************************************************************************************
## Overview

`Dornic_et_al_integration_method` provides a general class to integrate stochastic differential equations with demographic multiplicative noise.

It is an optimization and generalization of a basic code implemented by Ivan Dornic and Juan A. Bonachela based on the algorithm proposed in:

> “Integration of langevin equations with multiplicative noise and the viability of field theories for absorbing  phase transitions”; Ivan Dornic, Hugues Chaté, and Miguel A Munoz;  Physical review letters 94; 100601 (2005).

and improved following:

> “Simulation of spatial systems with demographic noise”; Haim Weissmann, Nadav M Shnerb, and David A Kessler; Physical Review E 98; 022131 (2018).

This project is licensed under the terms of the GNU General Public License license. If you use our code in an academic publication, please include us in the acknowledgements and cite the previous papers. 

***********************************************************************************************************************************

## Installation and compiling

Download or clone this repository using `git clone`. Integration routines are inside the file`Dornic_et_al_integration_method.h`, that should included in your C++ file as usual. Then compile the code,

```bash
g++ *.cpp -std=c++11 -other_flags -o example_main.exe
```
The Dornic integration method provides also three compilation compilation options: 

- `CONSTANT_CELL_NUMBER`: indicates whether the number of nodes during the integration is constant. It defaults to `true`. If you want to integrate in a complex network with changing number of nodes, you should set it to `false`.
- `USE_EULER`: if set to `true`, the deterministic part is integrated using the Euler algorithm instead of the default RK4. By default it is set to `false`. Notice that the bottleneck of the algorithm is in the Poisson number generation in conjuction with  low `dt`(see "Limitations" section below) hence for many cases changing the method has little impact on execution time and is not recommended.
- `RNG`: indicates the random engine used to generate the random numbers. It should be set to match your generator of choice. By default it is the Mersenne Twister `mt19937`.
- `APPROX_POISSON`: when the mean of the Poissonian distribution is  larger than this macro, a normal distribution is used instead of the Poisson. If set to zero, it has no effect and all values will be distributed as a Poissonian. This is done because the main bottleneck is Poissonian number generation. Its default (and recommended) value is zero. Please see the discussion at "Limitations" section below before tweaking with this parameter.

An example compilation changing these parameters would be

```bash
g++ example_main.cpp -std=c++11 -other_flags -DCONSTANT_CELL_NUMBER=false -DRNG=mt19937_64 -o example_main
```

***********************************************************************************************************************************

##  How to use it

### Setting the equation

The header file `Dornic_et_al_integration_method.h` contains all the core integration routines, therefore it shouldn't be changed by the user for most applications.  The algorithm already takes care of the multiplicative noise and the linear part, so only linear differential operators (like diffusion) and non-linear parts have to be integrated separately. The user can specify the right hand side of the equation inside `Dornic_et_al_integration_method.cpp`, in the method `non_linear_rhs`.  So, if the equation reads $\dot\phi(x,t) = \mathcal L[\phi] + \mathcal N[\phi] +\mathcal D[\nabla\phi]+\sqrt{\phi}\eta(t)$, the user needs to implement the part corresponding to $\mathcal N [\phi]+ \mathcal D[\nabla\phi ]$. 

As an example, the right hand side of the Reggeon field theory has been already implemented ($\mathcal N[\phi]=-b\phi^2$, $\mathcal D[\nabla\phi]=D\nabla^2\phi​$, see also `Dornic_method.pdf`), as

```c++
double Dornic::non_linear_rhs(const int inode, const vector<double> &field) const //[1]
{
    //Non-linear terms [2]
    const double quadratic = -quad_coeff * field[inode] * field[inode];

    //Diffusion 
    double diff_sum;
    int num_neigh, index_neigh;

    num_neigh = neighbors[inode].size(); //[3]
    diff_sum = 0.0;

    //DIFFUSION INTEGRATION
    for(int i=0; i<num_neigh; i++)
    {
        index_neigh = neighbors[inode][i]; //[3]
        diff_sum += field[index_neigh];
    }

    diff_sum = D*(diff_sum - num_neigh*field[inode]);

    return diff_sum + quadratic; //[4]
}
```

Comments to this code are the following

1. The first line is the function specification, which should not be changed. As you can see, it belongs to the `Dornic` class (so you have access to every class private variable -be careful!) and receives an `inode` and a `field`. In general, you should consider `field[inode]` as the term $\phi(x,t)$ in your equation...
2. ...so, for example, the quadratic term $-b \phi^2$ of the Reggeon field theory translates to `-quad_coeff * field[inode]*field[inode]`. The definition of `quad_coeff` will be seen later.
3. The variable `neighbors` is defined inside the `Dornic` class. `neighbors[inode][j]` contains the index of the `j`-th neighbor of the node `inode`. Therefore, a network, the degree of `inode` can be go by `neighbors[inode].size()`.
4. The function must return a double, which is the result of the non-linear and differential operators part of your right hand side equation. 

The user-defined variables, as `quad_coeff` and `D`, can be set in the function `set_non_linear_coefficients`, in the same file. In order to do that, they have to be defined as global variables, and receive the information from the `f_parameters` argument in the function, as

```c++
//User-defined variables [1]
double quad_coeff;
double D;

void Dornic::set_non_linear_coefficients(const vector <double> &f_parameters)
{
    quad_coeff = f_parameters[2]; //[2]
    D = f_parameters[3] / (dx * dx); //[3]
}
```

The comments to this code are:

1. The variables are defined as global so they can be also used in the function `non_linear_rhs`. They are _not_ stored inside the Dornic class at any means.
2. The vector `f_parameters` is passed by reference to the Dornic class when this is constructed. **The first two parameters must be always the linear and noise coefficients, respectively**. For this reason, the first 'custom' parameter free for the user is `f_parameters[2]`.
3. The variable `dx` is defined inside the Dornic class. For a list of all variables inside the class, you always can refer directly to the `Dornic_et_al_integration_method.h` file, at the list of private variables. 

### Integrating with Dornic _et al_ package

In order to show the proper way of using it, we implemented the `example_main.cpp`. These codes integrate  the Langevin equation of the Reggeon field theory described in `Dornic_method.pdf` in a 2D lattice.  

The first is to declare the Dornic package, indicating the timesteps for both space and time, as well as the number of nodes in the system. The coefficients can be set in the constructor or later. 

```c++
//Parameters: linear and noise go always first!!
vector<double> f_parameters = {linear_coefficient, noise_coefficient, quadratic_coefficient, diffusion_coefficient};

Dornic dornic(dt,dx,number_cells);
dornic.set_coefficients(f_parameters);

//Can be done in one step as
//Dornic dornic(dt,dx,number_cells,f_parameters);
```

> **A very important remark is that the first two parameters must be the linear and noise coefficients, in that order. After those, the user can write parameters in the desired order.** 

Then,  the next step is to create an adequate lattice. The `Dornic` class has three built-in method for this task,  `set_1D_lattice`, `set_2D_lattice` and `custom_network`. The first two are self-explanatory, creating a regular lattice in one or two dimensions, respectively. The last one takes any `vector<vector<int>>` as argument, to create any network geometry you wish. Topology is independent from dynamics, so this function should be called just once.

For initial conditions, there are three functions for the most common use-cases: 

- `random_initial_cond(min,max)` gives a random value between `min` and `max` to each cell. If not specified, `min=0` and `max=1` by default.
- `homogeneous_initial_cond(value)` sets everybody to the same density value. If this is not specified, takes by default `value=1`.
- `single_seed(inode,value)` sets the density of the specified `inode` to `value`. All others are zero. If `value` is not set, it will be `value=1` by default.

However, if you wish to use any other initial conditions, this is also possible by changing manually the values of the density cells. The density of cell `inode` can be accessed (and modified) by calling `dornic[i]`. For example, making only odd cells to have activity 1 would go as

```c++
Dornic dornic(dt,dx,number_cells, f_parameters);

dornic.set_1D_lattice();

for (i=0; i < N; i++) dornic[i] = i%2 ? 1.0 : 0.0;
```

The average activity can be retrieved (but not modified) by calling `dornic.density()`.  Note that manually changing the cells will not update the total density until an integration step is performed. In general, tweaking densities between integration steps is not recommended.

Finally, to integrate, we just call `dornic.integrate(rng)`, which performs _one_ integration step. For the case where the coefficients are itself functions of time, `dornic.integrate(rng, f_parameters)` allows to change the parameters each integration steps.  The variable `rng` is our random number generator, which is recommended to be of type `RNG` to coincide with the one specified in compilation options (see above). For example, the following code makes a simple integration during a time interval

```c++
RNG gen(seed); //Initialize the RNG using the pre-compiled type

//Dornic initialization...

ft_file.open(filepath);
for (t=0; t<=tmax; t+=dt)
{
    dornic.integrate(gen); //Integration step
    ft_file << t << " " << dornic.density() <<endl;
}
ft_file.close();
```

For a complete, functional code using Dornic integration, please refer to the `example_main.cpp` file, where two examples are implemented: 

1. Simple integration, writing total density vs time.
2. Simple phase diagram estimation of the contact process, running the program multiple times and computing average and variance of total density.

## Limitations

The code presented here still has some limitations. Let us speak a bit about an important one and how we tackled it.

The bottleneck of the code is the Poisson number generation. The linear parameter controls the mean of the distribution, $\mu$. When the linear parameter goes to 0, the mean of the Poissonian distribution grows as `noise_coeff * site_density / dt`, which can be large if the local site density is high and we want low `dt` values. The expected time to generate a Poissonian random variable grows with its mean. The C++ STL treats efficiently this case, but, since the mean of the distribution has to be re-computed each iteration (because it depends on local cell density) the generation is still very slow, due to many calls to `log`, `sqrt`, and `lgamma` functions. In order to alleviate this issue when needed, we included the `APPROX_POISSON` macro. When it is set to 0 (default, recommended), the Dornic algorithm uses the real Poisson distribution, no matter how much its mean is. However, when `APPROX_POISSON > 0`, if the mean is larger than the value indicated by `APPROX_POISSON`, the Poissonian distribution is substituted by a normal distribution, since in the limit of large mean, the Poisson behaves as a normal distribution with mean $\mu$ and variance $\mu$. The relative error of this approximation is bounded by $\max |G-P|/P =1/(3\sqrt\mu)​$, as shown by [Curtis](https://aapt.scitation.org/doi/10.1119/1.10057). 

In general, the best way to avoid this limitation is letting `dt` to be as higher as possible (that's why the integration scheme is a fixed time Runge-Kutta, which allows fairly large timesteps). You always can inspect the average value that the Poisson distribution is receiving with `dornic.avg_poisson_mean()`. If these values are very high and you _really_ need more speed, you should try to set the `APPROX_POISSON` macro to cut some values. Be warned that cutting low values could lead to inaccuracies in your results

Another usability limitation is for systems of equations. The code is though (by now) to integrate a single field equation subject to demographic noise. Multiple field could include different noise sources, increasing a lot the number of possible cases -so they have not been included here. The case for non-autonomous systems, which depend explicit on time is not yet supported, but we expect they will be. 

Finally, _efficient_ mean-field computations are still not possible with this program. Of course, you always can simulate a fully-connected topology, but this contains a lot of calls to neighbors, hence resulting in very slow code. Mean-field computations can be made better by just letting each particle interact with the mean field. This will be supported in the future.

***********************************************************************************************************************************