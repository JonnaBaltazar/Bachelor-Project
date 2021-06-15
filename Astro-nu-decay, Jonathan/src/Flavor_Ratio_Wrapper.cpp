#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <stdio.h>
#include <vector>
#include <math.h>

#include "Ratio.h"
#include "Parameters.h"
#include "Figures.h"
#include "Redshift.h"
#include "Energy_Integration.h"
#include "IceCube.h"
#include "SM.h"
#include "Functions.h"

/* Written by Mauricio Bustamane & Jonathan Baltazar 
email = mbustamante@nbi.ku.dk
email = jonatuk@gmail.com 

Last modified: 11/6-2021
*/

/* Add all the includes you need here */
Ratio flavor_wrapper(double flavor_source[3], double mnu1, double gamma, \
    double g21, double g31, double g32, double s12sq, double s23sq, double s13sq, double deltacp, int ykbh)
{	
    decay_params dp = dp_bm; 
    mix_params mp = mp_bm; 

    dp.m1 = mnu1;
    dp.gamma = gamma;
    dp.g21 = g21;
    dp.gp21 = g21;
    dp.g31 = g31;
    dp.gp31 = g31;
    dp.g32 = g32;
    dp.gp32 = g32;

    mp.s12sq = s12sq;
    mp.s23sq = s23sq; 
    mp.s13sq = s13sq;
    mp.delta = deltacp; 

    Ratio flavor_earth; 
    
  //  std::cout << "Wrapper: "<< ykbh << std::endl;
    
    double Emin=6*1e4, Emax=3*1e6;

    flavor_earth = integrate_flavor(Emin, Emax, flavor_source, dp, ykbh); 

    std::cout << flavor_earth.Fe << " "  << flavor_earth.Fmu << " "  << flavor_earth.Ftau << std::endl;
    
    return flavor_earth; 
}

/* Defining all the parameters that the function must take into account: 
mnu1_npts=4 other =1.
num_blocks=2, index_block=0  

*/
void generate_flavor_ratios_earth_block( int index_block, int num_blocks, \
    double mnu1_min, double mnu1_max, int mnu1_npts, \
    double gamma_min, double gamma_max, int gamma_npts, \
    double log10g21_min, double log10g21_max, int log10g21_npts, \
    double log10g31_min, double log10g31_max, int log10g31_npts, \
    double log10g32_min, double log10g32_max, int log10g32_npts, \
    double s12sq_min, double s12sq_max, int s12sq_npts, \
    double s23sq_min, double s23sq_max, int s23sq_npts, \
    double s13sq_fixed, double deltacp_fixed, double F_e_s, double F_mu_s,\
    char* output_path, int verbose, int ykbh)
{ 
    int num_points_total, num_points_in_block;
    int index_start, index_end;

    int mnu1_index;
    int gamma_index;
    int log10g21_index, log10g31_index, log10g32_index;
    int s12sq_index, s23sq_index;
    int k1, k2, k3, k4, k5, k6;

    double mnu1_step;
    double gamma_step;
    double log10g21_step, log10g31_step, log10g32_step;
    double s12sq_step, s23sq_step;

    double mnu1;
    double gamma;
    double log10g21, log10g31, log10g32;
    double g21, g31, g32;


    char filename[500], filename_full[500];
    std::ofstream output_file;  
 
    /* Defining the steplength for each of the 7 parameters that we want to vary: 
    step_length = (para_max - para_min) / (npt - 1) */
    mnu1_step = (mnu1_max-mnu1_min) / ((double)mnu1_npts-1.0);
    gamma_step = (gamma_max-gamma_min)/((double)gamma_npts-1.0);
    log10g21_step = (log10g21_max-log10g21_min)/((double)log10g21_npts-1.0);
    log10g31_step = (log10g31_max-log10g31_min)/((double)log10g31_npts-1.0);
    log10g32_step = (log10g32_max-log10g32_min)/((double)log10g32_npts-1.0);
    s12sq_step = (s12sq_max-s12sq_min)/((double)s12sq_npts-1.0);
    s23sq_step = (s23sq_max-s23sq_min)/((double)s23sq_npts-1.0);

    // Save the values of the parameters varied in external files, once only
    if (index_block == 0)
    {
		const char* Header_1 = "Parameter Values:";
		const char* Header_2 = "Collumn1=para_min, Collumn2=para_max, Collumn3=para_npts, Collumn4=para_step:";
		const char* Header_3 = "mnu1 --> gamma --> log10g21 --> log10g31 --> log10g32 --> s12sq --> s23sq:";

    	FILE *data = fopen("data/Parameter_Values.txt", "w");
    	setbuf(data, NULL);
    	fprintf(data, "%s\n%s\n%s\n\n", Header_1, Header_2, Header_3);
		fprintf(data, "%g %g %i %g\n", mnu1_min, mnu1_max, mnu1_npts, mnu1_step);
		fprintf(data, "%g %g %i %g\n", gamma_min, gamma_max, gamma_npts, gamma_step);
		fprintf(data, "%g %g %i %g\n", log10g21_min, log10g21_max, log10g21_npts, log10g21_step);
		fprintf(data, "%g %g %i %g\n", log10g31_min, log10g31_max, log10g31_npts, log10g31_step);
		fprintf(data, "%g %g %i %g\n", log10g32_min, log10g32_max, log10g32_npts, log10g32_step);
		fprintf(data, "%g %g %i %g\n", s12sq_min, s12sq_max, s12sq_npts, s12sq_step);
		fprintf(data, "%g %g %i %g\n", s23sq_min, s23sq_max, s23sq_npts, s23sq_step);
		fclose(data);
    }
    
    /* Compute the start and end indices of the parameter combinations to run
       for this particular block */

    // Total number of points: 
    num_points_total = mnu1_npts*gamma_npts \
        *log10g21_npts*log10g31_npts*log10g32_npts \
        *s12sq_npts*s23sq_npts;
    
    // Number of points in each block: 
    num_points_in_block = num_points_total/num_blocks;
    
    // We start at the number of points we have in each block, times by the block number we want to run: 
    index_start = index_block*num_points_in_block;

    // We stop before 1: 
    index_end = index_block*num_points_in_block+(num_points_in_block-1);
    
    if (index_block == num_blocks-1)
    {   
        // The % operater finds the modulus / remainder between two numbers: 
        index_end += num_points_total%num_blocks;
    }

    if (verbose)
    {
        std::cout << "Computing flavor ratios for index_block = " << index_block << std::endl;
    }
    if (ykbh)
    {
        std::cout << "Using the YKBH distribution!" << std::endl; 
    }

    /* The k numbers are used below to compute the values of the parameters in
       every iteration of the loop
       k0 = 1; */ 

    k1 = s23sq_npts;
    k2 = k1*s12sq_npts;
    k3 = k2*log10g32_npts;
    k4 = k3*log10g31_npts;
    k5 = k4*log10g21_npts;
    k6 = k5*gamma_npts;

    // Generate output file name and open it
    output_file << std::scientific << std::setprecision(5);
    sprintf(filename, "nu_flux_ratios_block_%d.dat", index_block);
    strcat(strcat(strcpy(filename_full, output_path), "/"), filename);
    output_file.open(filename_full, std::ios::out | std::ios::trunc);

    s13sq = s13sq_fixed; 
    delta = deltacp_fixed; 
    
    I2_precalc();

    // Generate and save the neutrino flux for each of the parameter 
    // combinations in this block
    for (int i=index_start; i<=index_end; i++)
    {
        if (verbose)
        {
            // std::cout << i-index_start \
            //     << "/" << num_points_in_block-1 << std::endl;
            std::cout << i-index_start << "/" << index_end-index_start << std::endl;
        }

        // std::cout << filename_full << std::endl;

        // Find the parameter values for this value of the index i
        // Note: the order in which parameters are varied is
        // mnu1 --> gamma --> log10g21 --> log10g31 --> log10g32 --> s12sq --> s23sq
        
        /* The conditional operator evaluates an expression, returning one value if that expression evaluates to true, and a different one if the expression evaluates as false. Its syntax is:

           condition ? result1 : result2

           If condition is true, the entire expression evaluates to result1, and otherwise to result2. */

        mnu1_index = (mnu1_npts > 1) ? (i/k6)%mnu1_npts : 0;        

        // If mnu1_npts > 1 the mnu1_index will evaluate as the difference between (i/k6) and mnu1_npts. 
        // Same for the rest of the indexes: 

        gamma_index = (gamma_npts > 1) ? (i/k5)%gamma_npts : 0;        
        log10g21_index = (log10g21_npts > 1) ? (i/k4)%log10g21_npts : 0;        
        log10g31_index = (log10g31_npts > 1) ? (i/k3)%log10g31_npts : 0;        
        log10g32_index = (log10g32_npts > 1) ? (i/k2)%log10g32_npts : 0;        
        s12sq_index = (s12sq_npts > 1) ? (i/k1)%s12sq_npts : 0;        
        s23sq_index = (s23sq_npts > 1) ? i%s23sq_npts : 0;        

        // If npts > 1 mnu1 will be evaluated as mnu1_min + the length of the steps and the index of mnu1. 
        mnu1 = (mnu1_npts > 1) ? mnu1_min + mnu1_step*mnu1_index : mnu1_min;     
        
        // This is the same as for the rest of the parameters: 
        gamma \
            = (gamma_npts > 1) ? gamma_min + gamma_step*gamma_index : gamma_min;     
        log10g21 = (log10g21_npts > 1) ? log10g21_min \
            + log10g21_step*log10g21_index : log10g21_min;     
        log10g31 = (log10g31_npts > 1) ? log10g31_min \
            + log10g31_step*log10g31_index : log10g31_min;     
        log10g32 = (log10g32_npts > 1) ? log10g32_min \
            + log10g32_step*log10g32_index : log10g32_min;     
        s12sq = \
            (s12sq_npts > 1) ? s12sq_min + s12sq_step*s12sq_index : s12sq_min;     
        s23sq = \
            (s23sq_npts > 1) ? s23sq_min + s23sq_step*s23sq_index : s23sq_min;  

        // Finally we can revert the original angles: 
        g21 = pow(10.0, log10g21);
        g31 = pow(10.0, log10g31);
        g32 = pow(10.0, log10g32);
        
        /* Compute the flavor ratios (modify to your taste).  This function
         * returns the result in the structure flavor_ratios */

        double F_tau_s;
        F_tau_s = 1 - F_e_s - F_mu_s;
        double init_ratio[3] = {F_e_s, F_mu_s, F_tau_s};

        Recalc_Parameters();

        Ratio final_ratio; 

        final_ratio = flavor_wrapper(init_ratio, mnu1, gamma, g21, g31, g32, s12sq, s23sq, s13sq_fixed, deltacp_fixed, ykbh);

        /* Save result: first column is fe, second column is fmu */
        output_file << final_ratio.Fe << " " << final_ratio.Fmu << std::endl;
    }

    output_file.close();
}

/*  The main can be changed to this by rewriting 
    void int_main to int main
    and removing // before return 0; 

    The other main function, from main.cpp, just needs to be commented out.
*/
int main(int argc, char **argv)
{	

    int index_block, num_blocks; 

    int log10g21_npts, log10g31_npts, log10g32_npts;
    int gamma_npts;
    int mnu1_npts;
    int s12sq_npts, s23sq_npts;

    unsigned int verbose;
    int ykbh;

    double log10g21_min, log10g21_max;
    double log10g31_min, log10g31_max;
    double log10g32_min, log10g32_max;
    double gamma_min, gamma_max;
    double mnu1_min, mnu1_max;
    double s12sq_min, s12sq_max;
    double s23sq_min, s23sq_max;
    double F_e_s, F_mu_s;

    double s13sq, deltacp;

    char* output_path;

    // Read in the arguments passed in the command line
    // argv[0] is the file name
    index_block = atoi(argv[1]);
    num_blocks = atoi(argv[2]);
    mnu1_min = atof(argv[3]);
    mnu1_max = atof(argv[4]);
    mnu1_npts = atoi(argv[5]);
    gamma_min = atof(argv[6]);
    gamma_max = atof(argv[7]);
    gamma_npts = atoi(argv[8]);
    log10g21_min = atof(argv[9]);
    log10g21_max = atof(argv[10]);
    log10g21_npts = atoi(argv[11]);
    log10g31_min = atof(argv[12]);
    log10g31_max = atof(argv[13]);
    log10g31_npts = atoi(argv[14]);
    log10g32_min = atof(argv[15]);
    log10g32_max = atof(argv[16]);
    log10g32_npts = atoi(argv[17]);
    s12sq_min = atof(argv[18]);
    s12sq_max = atof(argv[19]);
    s12sq_npts = atof(argv[20]);
    s23sq_min = atof(argv[21]);
    s23sq_max = atof(argv[22]);
    s23sq_npts = atof(argv[23]);
    s13sq = atof(argv[24]);
    deltacp = atof(argv[25]);
    F_e_s = atof(argv[26]);
    F_mu_s = atof(argv[27]);
    output_path = argv[28];
    verbose = atoi(argv[29]);
    ykbh = atoi(argv[30]);

    generate_flavor_ratios_earth_block( \
        index_block, num_blocks, 
        mnu1_min, mnu1_max, mnu1_npts,
        gamma_min, gamma_max, gamma_npts,
        log10g21_min, log10g21_max, log10g21_npts,
        log10g31_min, log10g31_max, log10g31_npts,
        log10g32_min, log10g32_max, log10g32_npts,
        s12sq_min, s12sq_max, s12sq_npts,
        s23sq_min, s23sq_max, s23sq_npts,
        s13sq, deltacp,
        F_e_s, F_mu_s,
        output_path, verbose, ykbh);

    return 0;
}

