
// Written by Jonathan Baltazar (email: Jonatuk@gmail.com). 
// Last modified: 11/6-2021

#include <stdio.h>
#include <cmath>
#include <array>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#include "Energy_Integration.h"
#include "SM.h"
#include "Parameters.h"
#include "Depletion.h"
#include "Regeneration.h"
#include "Analytic.h"
#include "Redshift.h"
#include "Ratio.h"

// Three different functions to integrate, for each of the different flavors: 
double integrate_electron(double Ef, double ratio[3], decay_params dp, int ykbh)
{
	// For visible-decay scenario, depletion and regeneration are set to true: 
	bool include_dep = true;
	bool include_reg = true;

	// std::cout << "Electron: "<< ykbh << std::endl;
	// ykbh = 0; 
	// std::cout << "Electron after declaration: "<< ykbh << std::endl;
	
	double result;

	dp.Ef = Ef;

	Ratio fr_earth; 

	fr_earth = ratio_earth(ratio, dp, include_dep, include_reg, ykbh);
	
	result = fr_earth.Fe; 
	return result; 
}

double integrate_muon(double Ef, double ratio[3], decay_params dp, int ykbh)
{
	bool include_dep = true;
	bool include_reg = true;

	double result;
	dp.Ef = Ef;
	
	// std::cout << "Muon: "<< ykbh << std::endl;
	// ykbh = 0; 
	// std::cout << "Muon after declaration: "<< ykbh << std::endl;
	
	Ratio fr_earth; 
	fr_earth = ratio_earth(ratio, dp, include_dep, include_reg, ykbh);

	result = fr_earth.Fmu; 
	return result; 
}

// We define an integrand: 
double integrand(double Ef, void *params)
{
	// Using the structure from Energy_Integration.h define the parameters: 
    integrand_energy *p = (integrand_energy *)params;
   
	// Decay parameters. OBS: Missing are bool x_is_redshift and struct nunubar:
	double m1 = p->m1;
	double g21 = p->g21;
	double g31 = p->g31;
	double g32 = p->g32;
	double gp21 = p->gp21;
	double gp31 = p->gp31;
	double gp32 = p->gp32;
	double gamma = p->gamma;
	double x = p->x;


    double* ratio=p->ratio;
    int ykbh=p->ykbh;

	double (*flavor_function)(double Ef, double ratio[3], decay_params dp, int ykbh)=p->flavor_function;

	// All the parameters are initialized: 
	decay_params dp = dp_bm;
	dp.m1 = m1;
	dp.g21 = g21;
	dp.g31 = g31;
	dp.g32 = g32;
	dp.gp21 = gp21;
	dp.gp31 = gp31;
	dp.gp32 = gp32;
	dp.Ef = Ef;
	dp.gamma = gamma;
	dp.x = x;

    return flavor_function(Ef, ratio, dp, ykbh);
}

/* This function integrates the flavour function ratio_earth, over a given energy range and normalizes the results. This is the same way they 
   treat their data at IceCube and this makes sure that the simulated results are comparable with the actual ones. 
   Both the decay and mixing parameters can be varied, but it is set to visible decay by default - including both regeneration and depletion. 
   The accompanying function Integrate_Flavor can be found in the file Integrate_Flavor.cpp. */

Ratio integrate_flavor(double Emin, double Emax, double ratio[3], decay_params dp, int ykbh)
{	

    double result_electron;
    double result_muon;

    size_t neval;
    //int n = 50; 
	
    gsl_integration_romberg_workspace* w = gsl_integration_romberg_alloc(ROMB_N); // ROMBERG 
    // const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre; // Fixed Quadrature
    // gsl_integration_fixed_workspace* w = gsl_integration_fixed_alloc(T, n, Emin, Emax, 0.0, 0.0); 
    
    //gsl_integration_workspace* w = gsl_integration_workspace_alloc(ALLOC); // QAGS
    integrand_energy params;

    gsl_function F;
    F.function = &integrand;

	params.ratio[0] = ratio[0];
	params.ratio[1] = ratio[1];
	params.ratio[2] = ratio[2];

	params.m1 = dp.m1;
	params.g21 = dp.g21;
	params.g31 = dp.g31;
	params.g32 = dp.g32;
	params.gp21 = dp.gp21;
	params.gp31 = dp.gp31;
	params.gp31 = dp.gp32;
	params.Ef = dp.Ef;
	params.gamma = dp.gamma;
	params.x = dp.x;
	params.ykbh = ykbh; 

	double norm = 1 / (Emax-Emin);

	params.flavor_function = integrate_electron;
    F.params = &params;
    gsl_integration_romberg(   &F, Emin, Emax, EPSABS, EPSREL,
                               &result_electron, &neval, w);
    //gsl_integration_fixed(&F, &result_electron, w);
    // gsl_integration_qags(   &F, Emin, Emax, EPSABS, EPSREL, LIMIT, w,
    // 	                    &result_electron, &error_electron);    

    params.flavor_function = integrate_muon;
    F.params = &params;

    gsl_integration_romberg(   &F, Emin, Emax, EPSABS, EPSREL,
                               &result_muon, &neval, w);
    //gsl_integration_fixed(&F, &result_muon, w);
    // gsl_integration_qags(   &F, Emin, Emax, EPSABS, EPSREL, LIMIT, w,
    // 	                    &result_muon, &error_muon);

    
    gsl_integration_romberg_free(w); // Romberg 
    // gsl_integration_fixed_free(w); // Fixed Quadrature
    // gsl_integration_workspace_free(w); // QAGS


    Ratio result; 

    result.Fe = norm*result_electron; 
    result.Fmu = norm*result_muon; 
    result.Ftau = 1 - result.Fe - result.Fmu; 

    return result; 
}