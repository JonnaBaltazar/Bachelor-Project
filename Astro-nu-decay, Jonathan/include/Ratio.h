
// Written by Jonathan Baltazar (email: Jonatuk@gmail.com). 
// Last modified: 11/6-2021

#ifndef Ratio_Earth_H
#define Ratio_Earth_H

#include "Parameters.h"

// double Ratio_Earth(double fr_source;  struct mix_params, struct decay_params; bool include_dep, include_reg);

Ratio ratio_earth(double flavour_source[3], decay_params dp, bool include_dep, bool include_reg, int ykbh);

void generate_flavor_ratios_earth_block( int index_block, int num_blocks, \
    double mnu1_min, double mnu1_max, int mnu1_npts, \
    double gamma_min, double gamma_max, int gamma_npts, \
    double log10g21_min, double log10g21_max, int log10g21_npts, \
    double log10g31_min, double log10g31_max, int log10g31_npts, \
    double log10g32_min, double log10g32_max, int log10g32_npts, \
    double s12sq_min, double s12sq_max, int s12sq_npts, \
    double s23sq_min, double s23sq_max, int s23sq_npts, \
    double s13sq, double deltacp, double F_e_s, double F_mu_s,\
    char* output_path, int verbose);

#endif