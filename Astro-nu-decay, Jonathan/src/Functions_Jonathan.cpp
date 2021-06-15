
// Written by Jonathan Baltazar (email: Jonatuk@gmail.com). 
// Last modified: 11/6-2021

#include <stdio.h>
#include <cmath>
#include <array>
#include <iostream>

#include "Functions.h"
#include "SM.h"
#include "Energy_Integration.h"
#include "Depletion.h"
#include "Regeneration.h"
#include "IceCube.h"
#include "Analytic.h"
#include "Redshift.h"
#include "Ratio.h"

/* Calculate the flavor ratio at earth, given the flavour-ratio at the source. In this function decay parameters can be varied.
	It uses the function called ratio_earth, which is defined in the file Ratio_Earth.cpp with the accompanying header file Ratio.h. 
	It is normalized by default. */

void Flavor_Ratio_Earth()
{
	// Define the ratio at the source: 
	double fr_source[3] = {1./3.,2./3.,0};

	/* Choose wheather or not depletion and regeneration is going to be included. 
	Standard  = false, false 
	Invisible = true, false 
	Visible   = true, true
	*/

	bool include_dep = true;
	bool include_reg = true;

	// Decay parameters from Parameters.h and Parameters.cpp: 
	// decay_params {double m1, g21, g31, g32, gp21, gp31, gp32, Ef, gamma, x; bool x_is_redshift; nunubar nnb;}; 
	decay_params dp = dp_bm;


	// Turn on or off the couplings: 
	// dp.g21  = 0;
	// dp.gp21 = 0;
	// dp.g32  = 0;
	// dp.gp32 = 0;
	// dp.g31  = 0;
	// dp.gp31 = 0;

	// Choose delta function at z=1 or YKBH distribution: 
	int ykbh = 0;

	// Define the structure, from Ratio.h. The structure has three parameters: Fe, Fmu, Ftau 
	Ratio fr_earth; 

	// Preforming the calculation, normalized: 
	fr_earth = ratio_earth(fr_source, dp, include_reg, include_dep, ykbh); 

	printf("\nFlavor ratio at earth: \n");
	printf("1/3(1:2:0) -> (%.3g:%.3g:%.3g)\n", fr_earth.Fe, fr_earth.Fmu, fr_earth.Ftau);
}

/* This function calculates the flavour ratio as a function of the final energy, Ef, and stores the result in a data file called Flavor_Ef.txt,
   which can be found in the data folder.   
    */

void Flavor_Ef()
{
	printf("\nFlavor as a function of final Energy Ef\n");
	double fr_source[3] = {1./3.,2./3.,0};

	// The function calculates and stores the values for both visible decay and for the SM: 
	int ykbh = 0; 

	decay_params dp = dp_bm;

	Ratio fr_earth_in; // (include)
	Ratio fr_earth_ex; // (exclude)


	double Efmin, Efmax, Efscale;
	int n;

	// Number of iterations: 
	n = 1e2;

	// Initial and final Ef: 
	Efmin = 6*1e4;
	Efmax = 3*1e6;

	// Turn on or off the couplings: 
	// dp.g21  = 0;
	// dp.gp21 = 0;
	// dp.g31 = 0;
	// dp.gp31 = 0;
	// dp.g32 = 0;
	// dp.gp32 = 0;

	// A logarithmic scale for the energy:  
	Efscale = pow(Efmax / Efmin, 1. / n);

	// Creating file: 
	FILE *data = fopen("data/Flavor_Ef.txt", "w");
	setbuf(data, NULL);

	for (int i = 0; i <= n; i++)
	{
		printf("%g\n", 1.0 * i / n);
		dp.Ef = Efmin * pow(Efscale, i);
		
		fr_earth_in = ratio_earth(fr_source, dp,  true, true, ykbh); 
		fr_earth_ex = ratio_earth(fr_source, dp,  false, false, ykbh);

		// The data set has | Ef | F_e | F_mu | F_tau | F_e_SM | F_mu_SM | F_tau_SM |
		fprintf(data, "%g %g %g %g %g %g %g\n", dp.Ef, fr_earth_in.Fe, fr_earth_in.Fmu, fr_earth_in.Ftau, fr_earth_ex.Fe, fr_earth_ex.Fmu, fr_earth_ex.Ftau);
	}
	fclose(data);
} 

// Final flavor ratio as a function of energy, for different coupling constants:
void Flavor_Coupling_Ef()
{
	printf("Flavor_Coupling_Ef\n");

	double Efmin, Efmax, Efscale;
	int n;

	const int n_logg21 = 4;
	double logg21[n_logg21] = {-10, -6, -5.5, -2.5}; // the different spectral indices to use

	int ykbh = 0;

	double fr_source[3]={1./3., 2./3., 0};

	decay_params dp = dp_bm;

	// Turn on or off the couplings: 
	dp.g21 = 0;
	dp.gp21 = 0;
	dp.g31 = 0;
	dp.gp31 = 0;
	dp.g32 = 0;
	dp.gp32 = 0;

	n = 1e2;
	Efmin = 1e3;
	Efmax = 1e9;
 
	Efscale = pow(Efmax / Efmin, 1. / n);

	FILE *data = fopen("data/Visible_coupling_g21.txt", "w");
	setbuf(data, NULL);
	
	fprintf(data, "%g %g %g %g %g %i\n", dp.g21, dp.g31, dp.g32, dp.x, dp.m1, n_logg21);

	for (int i = 0; i < n_logg21; i++)
		fprintf(data, "%g ", logg21[i]);


	fprintf(data, "\n");
	for (int i = 0; i <= n; i++)
	{
		printf("%g\n", 1.0 * i / n);

		dp.Ef = Efmin * pow(Efscale, i);

		for (int j = 0; j < n_logg21; j++)
		{	
			Ratio standard; 
			standard =  ratio_earth(fr_source, dp,  false, false, ykbh);
			if (j == 0) fprintf(data, "%g %g %g ", dp.Ef, standard.Fe, standard.Fmu);

			// Turn on or off the couplings - opposite from above: 
			dp.g21  = pow(10, logg21[j]);
			dp.gp21 = pow(10, logg21[j]);
			// dp.g31  = pow(10, logg21[j]);
			// dp.gp31 = pow(10, logg21[j]);
			// dp.g32  = pow(10, logg21[j]);
			// dp.gp32 = pow(10, logg21[j]);

			Ratio p; 
			p = ratio_earth(fr_source, dp,  true, true, ykbh);
			fprintf(data, "%g %g ", p.Fe, p.Fmu);
		} 
		fprintf(data, "\n");
	}
	fclose(data);
}

// Final flavor ratio as a function of gamma:
void Flavor_Gamma()
{
	printf("\nFlavor as a function of spectral index, gamma\n");
	double fr_source[3] = {1./3.,2./3.,0};

	int ykbh = 0; 

	decay_params dp = dp_bm;

	Ratio fr_earth_in;
	Ratio fr_earth_ex;

	double gamma_min, gamma_max, gamma_scale;
	int n;

	// Turn on or off the couplings: 
	// dp.g21 = 0;
	// dp.gp21 = 0;
	// dp.g31 = 0;
	// dp.gp31 = 0;
	// dp.g32 = 0;
	// dp.gp32 = 0;

	n = 1e2;
	gamma_min = 1;
	gamma_max = 5;

	gamma_scale = pow(gamma_max / gamma_min, 1. / n);

	FILE *data = fopen("data/gamma_All.txt", "w");
	setbuf(data, NULL);
	for (int i = 0; i <= n; i++)
	{
		printf("%g\n", 1.0 * i / n);

		dp.gamma = gamma_min * pow(gamma_scale, i);
		
		dp.Ef = 1e6;
		fr_earth_in = ratio_earth(fr_source, dp,  true, true, ykbh); 
		fr_earth_ex = ratio_earth(fr_source, dp,  false, false, ykbh);
		
		fprintf(data, "%g %g %g %g %g %g %g\n", dp.gamma, fr_earth_in.Fe, fr_earth_in.Fmu, fr_earth_in.Ftau, fr_earth_ex.Fe, fr_earth_ex.Fmu, fr_earth_ex.Ftau);
	}
	fclose(data);
} 

// Final flavor ratio as a function of m1:
void Flavor_M1()
{
	printf("\nFlavor as a function of spectral index, mass m1\n");
	double fr_source[3] = {1./3.,2./3.,0};

	int ykbh = 0; 


	decay_params dp = dp_bm;

	Ratio fr_earth_in;
	Ratio fr_earth_ex;
	double m1_min, m1_max, m1_scale;
	int n;

	// Turn on or off the couplings: 
	// dp.g21 = 0;
	// dp.gp21 = 0;
	// dp.g31 = 0;
	// dp.gp31 = 0;
	// dp.g32 = 0;
	// dp.gp32 = 0;

	n = 1e2;
	m1_min = 1e-10;
	m1_max = 3.5;

	m1_scale = pow(m1_max / m1_min, 1. / n);

	FILE *data = fopen("data/Mass_All.txt", "w");
	setbuf(data, NULL);
	for (int i = 0; i <= n; i++)
	{	
		printf("%g\n", 1.0 * i / n);

		dp.m1 = m1_min * pow(m1_scale, i);
		
		dp.Ef = 1e6;
		fr_earth_in = ratio_earth(fr_source, dp,  true, true, ykbh); 
		fr_earth_ex = ratio_earth(fr_source, dp,  false, false, ykbh);

		fprintf(data, "%g %g %g %g %g %g %g\n", dp.m1, fr_earth_in.Fe, fr_earth_in.Fmu, fr_earth_in.Ftau, fr_earth_ex.Fe, fr_earth_ex.Fmu, fr_earth_ex.Ftau);
	}
	fclose(data);
}

/* Visible decay probabillity as a function of final energy for different redshift evolutions:
	This function closely resembles the Visible_R function in Figures.cpp - but has HERMES removed. */

void Visible_Redshift()
{
	printf("Visible Redshift\n");

	double Ef, Efmin, Efmax, Efscale, psm;

	const int n = 1e2;
	
	flavor a, b;
	double p_sm[n + 1], p_delta5[n + 1], p_delta1[n + 1], p_delta15[n + 1], p_ykbh[n + 1];
	
	Efmin = 1e3;
	Efmax = 1e9;
	Efscale = pow(Efmax / Efmin, 1. / n);

	// nu_mu to nu_e transition: 
	a = m;
	b = e;

	# pragma omp parallel for schedule(dynamic)
	for (int i = 0; i <= n; i++)
	{

		printf("%g\n", 1.0 * i / n);
		
		decay_params dp = dp_bm;

		dp.Ef = Efmin * pow(Efscale, i);

		// Delta functions redshift evolutions:
		dp.x = 0.5;
		psm = PSM(a, b, dp);
		p_delta5[i] = psm + Pdep(a, b, dp) + Preg(a, b, dp);
		dp.x = 1;
		p_sm[i] = PSM(a, b, dp);
		p_delta1[i] = p_sm[i] + Pdep(a, b, dp) + Preg(a, b, dp);
		dp.x = 1.5;
		psm = PSM(a, b, dp);
		p_delta15[i] = psm + Pdep(a, b, dp) + Preg(a, b, dp);

		// YKBH red shift evolution:
		p_ykbh[i] = RSE_integral(a, b, dp, true);
	}

	// Write to file
	FILE *data = fopen("data/Visible_Redshift.txt", "w");
	fprintf(data, "%i %i %g %g %g\n", a, b, dp_bm.g21, dp_bm.m1, dp_bm.gamma);
	for (int i = 0; i <= n; i++)
	{
		Ef = Efmin * pow(Efscale, i);
		fprintf(data, "%g %g %g %g %g %g\n", Ef, p_sm[i], p_delta5[i], p_delta1[i], p_delta15[i], p_ykbh[i]);
	}
	fclose(data);
}