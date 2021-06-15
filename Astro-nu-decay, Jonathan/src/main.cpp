// // #include <omp.h>

// #include "Parameters.h"
// #include "Figures.h"
// #include "Redshift.h"

// #include "IceCube.h"
// #include "SM.h"

// #include "Functions.h"

// // The main can either be set to this, or the main from Flavor_Ratio_Wrapper.cpp. 
// int main()
// {
// 	setbuf(stdout, NULL);

// 	//omp_set_num_threads(40);

// 	// Initializing functions
// 	Recalc_Parameters();
// 	I2_precalc();

// 	// Each command calculates the data for one figure, see src/Figures.cpp
// 	// Some take <1 second, some take O(10) minutes on one core, some take O(10) hours on 40 cores

// 	// Flavor_Ratio(); // quick
// 	// Invisible(); // quick
// 	// Visible_g();
// 	// Visible_m1();
// 	// Visible_gamma();
// 	// Visible_R(); // slow
// 	// Visible_f();
// 	// Rtc();
// 	// Analytic_Validate(); // quick
// 	// IC_gamma();
// 	// IC_gamma_2D(); // slow
// 	// Visible_SPS();
// 	// Visible_nnb();
// 	// Chisq(); // slowww

// 	// Functions from Functions_Jonathan.cpp: 
// 	// Modified by: Jonathan Baltazar (email=Jonatuk@gmail.com)
// 	// Last modified: 11/6-2021

// 	// Flavor_Ratio_Earth();
// 	Flavor_Ef(); 
// 	// Flavor_Coupling_Ef();
// 	// Flavor_Gamma(); 
// 	// Flavor_M1(); 
// 	// Visible_Redshift();
	
// 	return 0;
// }
