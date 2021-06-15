
// Written by Jonathan Baltazar (email: Jonatuk@gmail.com). 
// Last modified: 11/6-2021


#include <stdio.h>
#include <cmath>
#include <array>
#include <iostream>

#include "SM.h"
#include "Parameters.h"
#include "Depletion.h"
#include "Regeneration.h"
#include "Analytic.h"
#include "Redshift.h"

Ratio ratio_earth(double flavor_source[3], decay_params dp, bool include_dep, bool include_reg, int ykbh)
{
	Ratio fr_earth_norm;
	
	if (ykbh==0)
	{
		double norm; 
		Ratio fr_earth;

		// Parameters:
		double prob[3][3];

		double flavor_earth[3] = {0, 0, 0};

		for (int a = 0; a < 3; a++)
	    {
	    	for (int b = 0; b < 3; b++)
	    	{
	    		prob[a][b] = PSM(flavor(a), flavor(b), dp); // SM
	    	
	    		if (include_dep)
	            prob[a][b] += Pdep(flavor(a), flavor(b), dp); // Depletion

	        	if (include_reg)
	            prob[a][b] += Preg(flavor(a), flavor(b), dp); // Regeneration
	    	}
	    }

	    for (int a=0; a < 3; a++)
		{
			for (int b=0; b < 3; b++)
				flavor_earth[a] += flavor_source[b] * prob[b][a];
		}

		fr_earth.Fe = flavor_earth[0];
		fr_earth.Fmu = flavor_earth[1];
		fr_earth.Ftau = flavor_earth[2];

		norm = flavor_earth[0] + flavor_earth[1] + flavor_earth[2];
		fr_earth_norm.Fe = fr_earth.Fe/norm;
		fr_earth_norm.Fmu = fr_earth.Fmu/norm;
		fr_earth_norm.Ftau = fr_earth.Ftau/norm;
	}
	
	if (ykbh==1)
	{
		double norm_ykbh; 
		Ratio fr_earth_ykbh;

		// Parameters:
		double prob_ykbh[3][3];

		double flavor_earth_ykbh[3] = {0, 0, 0};

		for (int a = 0; a < 3; a++)
    	{
    		for (int b = 0; b < 3; b++)
    		{
    			prob_ykbh[a][b] = RSE_integral(flavor(a), flavor(b), dp, true);
    		}
    	}

    	for (int a=0; a < 3; a++)
		{
			for (int b=0; b < 3; b++)
				flavor_earth_ykbh[a] += flavor_source[b] * prob_ykbh[b][a];
		}

		fr_earth_ykbh.Fe   = flavor_earth_ykbh[0];
		fr_earth_ykbh.Fmu  = flavor_earth_ykbh[1];
		fr_earth_ykbh.Ftau = flavor_earth_ykbh[2];

		norm_ykbh = flavor_earth_ykbh[0] + flavor_earth_ykbh[1] + flavor_earth_ykbh[2];
		fr_earth_norm.Fe = fr_earth_ykbh.Fe/norm_ykbh;
		fr_earth_norm.Fmu = fr_earth_ykbh.Fmu/norm_ykbh;
		fr_earth_norm.Ftau = fr_earth_ykbh.Ftau/norm_ykbh;
	}

	return fr_earth_norm;
}