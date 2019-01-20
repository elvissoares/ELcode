/*
    General Description:
    ====================
 
        This header file contains the most up-to-date physical and
        astronomical constants in SI units.
        
        Confered in http://physics.nist.gov/cuu/Constants/index.html 
        in 10/28/15
 
 -------------------------------------------------------------------*/

#ifndef _H_CPP_CONSTANTS
#define _H_CPP_CONSTANTS

#include <cmath>

#define Q(x) ((x)*(x))
#define C(x) ((x)*(x)*(x))

namespace PhysConstants
{

//Values related to pi 
 const long double	pi          = acos(-1);
 const long double	two_pi      = 2*pi;
 const long double	four_pi     = 4*pi;
 const long double	four_pi_o3  = four_pi/3.;
 const long double	pi_over_2   = pi/2.;

//Conversions for radians to degrees and degrees to radians
 const long double  	d2rad 	= pi/180.;                 // degrees to radianes
 const long double  	rad2d 	= 180/pi;                // radianes to degrees

//Conersion constants

 const double       pc2m 	= 3.08568025e+16;       // 1 pc in meters
 const double       in2m 	= 0.0254;               // inches to meter conversion
 const double       MeV2g 	= 1.78358209e-27;   	// MeV to g conversion
 const double       GeV2kg 	= 1.78358209e-27;   	// GeV to kg conversion
 const double       MeV2erg = 1.60217653e-6;    	// MeV to erg
 const double       GeV2J 	= 1.60217653e-10;    	// GeV to Joule

// Universal constants
 const double       c       = 2.99792458e+10;       // speed of light in cm/s
 const double       G       = 6.67408e-8;        	// Gravitational cte (cm³/g.s²)
 const double		mu_0	= 1.25664e-6;          	// vacuum permeability (N/A²)
 const double		eps0 	= 8.854187817e-12;   	// vacuum permittivity (F/m)

 const double		h       = 6.626070040e-27;        // Planck cte (erg.s)
 const double		hbar	= h/two_pi;
 
 const double		k_B     = 1.380658e-16;        // Boltzmann cte (erg/K)
 const double		N_A     = 6.0221367e23;         // Avogadro number (mol^-1)

// Electromagnetic constants
 const double		e_c     = 4.80320425e-10;       // Electron charge (statcoulombs)

 const double		mu_B	= 9.2740154e-24;     	// Bohr Magneton (J/T)
 const double		mu_N	= 5.0507866e-27;     	// Nuclear Magneton (J/T)

// Atomics constants
 const double		alpha_f	= 7.29735308e-3;        // Fine structure constant
 const double		R_infty	= 1.0973731534e7;     	// Rydberg cte (m^-1)
 const double		a0      = 0.529177249e-10;     	// Bohr radius (m)

//Time constants
 const short int	hr      = 3600;
 const long int     day    	= 24*hr;
 const double       J_yr  	= 365.25*day;
 const double       yr     	= 3.15581450e7;  //year in seconds
 const double       T_yr   	= 3.155692519e7;
 const double       G_yr    = 3.1556952e7;

//Astronomical length constants
 const double        	AU    	= 1.4959787066e13;	// 1 AU in cm
 const double        	pc     	= 206264.806*AU;	// 1 pc in meters
 const double        	ly     	= c*yr;		// 1 light-year in cm

//Solar constants
 const float         	M_Sun  	= 1.98855e33F;		// Sun mass (g)
 const float         	R_Sun  	= 6.95508e10F;		// Sun radius (cm)
 const float         	L_Sun   = 3.9e33F;		// Luminosity (erg/s)
 const float         	T_Sun  	= 5778;			// Surface temperature of the Sun (K)

//Earth constants
 const float		M_Earth = 5.9742e27F;       	// Earth mass (g)
 const float		R_Earth	= 6.3781e8;        	// Earth radius (cm)
 const float	    rho_Earth	= 5.515;           	// Earth medium density (g/cm³)
 const float		g	= 980.62;          	// gravitational aceleration (cm/s²)
 const float 		P_atm	= 101.325;          	// atmospheric pressure (kPa) 

// Electron constants
 const double		M_e	= 9.10938356e-28;		// Electron mass (g)
 const double		R_e	= 2.81794092e-13; 	// Classical electron radius (cm) R_e=sqrt(e²/mc²)
 const double		sigma_e	= 0.66524616e-24; 	//Thomson cross section (cm²)
 const double		mu_e	= 9.2847701e-24;  	// Electron magnetic moment (J/T)

// Proton constants
 const double		M_p	= 1.6726231e-24;  	// Proton mass (g)
 const double		mu_p	= 1.41060761e-26;  	// Proton magnetic moment (J/T)

// Neutron
 const double		M_n	= 1.6749286e-24;  	// Neutron mass (g)
 const double		mu_n	= 0.96623707e-26;  	// Neutron magnetic moment (J/T)


 const float 		Tcmb	= 2.725;             	// CMB temperature today (K)

//Unit Conversions
    const float         m          	= 1e2;
    const float         kg          = 1e3;
    const float         joule       = 1e7;
    const float         dyne        = 1e-5;
    const double        esu         = 3.335640952e-10;
    const double        statvolt    = 2.997924580e2;
    const float         gauss       = 1e-4;
    const float         angstrom    = 1e-10;
    const float         jansky      = 1e-26; 
    const double 	    amu 	    = 1.660538921e-24;

    const double 		a 	= 7.5646e-15; //Black-body radiation constant (erg/cm^3K^4)
}

#endif
