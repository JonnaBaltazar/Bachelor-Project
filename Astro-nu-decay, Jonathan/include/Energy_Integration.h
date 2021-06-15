
// Written by Jonathan Baltazar (email: Jonatuk@gmail.com). 
// Last modified: 11/6-2021

#ifndef INTEGRATION_H
#define INTEGRATION_H

/* Target errors, memory, limit number of iterations for the GSL integration. Try 10^-3, 5*10^-3 and 10^-2. */
#define EPSABS    1.e-3 /* Absolute error */
#define EPSREL    1.e-3  /* Relative error */
#define ALLOC     1000   /* Allocated memory */
#define LIMIT     1000   /* Limit number of iterations */
#define ROMB_N    10     /* Used for Rombert integration */

#include "Parameters.h"

// Defining the integration functions: 
Ratio integrate_flavor(double Emin, double Emax, double ratio[3], decay_params dp, int ykbh);

// Defining flavor functions: 
double integrate_electron(double Ef, double ratio[3], decay_params dp, int ykbh);
double integrate_muon(double Ef, double ratio[3], decay_params dp, int ykbh);

// A structure for the Integrand is created: 
struct integrand_energy
{	

	// The parameters from decay_params is defined. OBS: Missing bool x_is_redshift and struct nunubar: 
	double m1;
	double g21;
	double g31;
	double g32;
	double gp21;
	double gp31;
	double gp32;
	double Ef;
	double gamma;
	double x; 

	int ykbh;
	// An array for the initial flavor ratio is defined: 
	double ratio[3];

	// The flavor function, of the same form as the integrate_flavor functions from above is defined: 
	double (*flavor_function)(double, double*, decay_params, int);
};

// Defining the structure: 
typedef struct integrand_energy integrand_energy;

#endif