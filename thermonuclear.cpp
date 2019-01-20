#include "thermonuclear.hpp"

#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "cpp_constants_CGS.hpp"
#include "util/nr.hpp"
#include "util/roots.hpp"

// This routine sets up the system of equations for the nuclear network
// input is y[1:7] the abundance fractions
// output is dydt[1:7] the right hand side of the network equations
//
// Isotopes: he4, c12, o16, ne20, mg24, si28, ni56

void Nuclear_Network::Calculate_Derivatives(const double &rho,
 const double &T, const std::vector<double> &X, std::vector<double> &dXdt)
{
    std::vector<double> Y(7);
    const double A[7] = {4.0,12.0,16.0,20.0,24.0,28.0,56.0};
    
    for (unsigned int i = 0; i < 7; i++)
        Y[i] = X[i]/A[i];

	Calculate_Rates(rho,T,Y);

	// 4He reactions
	dXdt[he4] = A[he4] * (3.0 * Y[c12] * rate[rg3a]
				- 3.0 * Y[he4] * Y[he4] * Y[he4] * rate[r3a]
				- Y[c12] * Y[he4] * rate[rcag]
				+ Y[c12] * Y[c12] * rate[r1212]
				+ 0.5 * Y[c12] * Y[o16] * rate[r1216]
                + Y[o16] * rate[roga]
				+ Y[o16] * Y[o16] * rate[r1616]
				- Y[o16] * Y[he4] * rate[roag]
				+ Y[ne20] * rate[rnega]
				- Y[ne20] * Y[he4] * rate[rneag]
				+ Y[mg24] * rate[rmgga]
				- Y[mg24] * Y[he4] * rate[rmgag]
				+ Y[si28] * rate[rsiga]
				- 7.0 * rate[rsi2ni] * Y[he4]
				+ 7.0 * rate[rni2si] * Y[ni56]);

	// 12C reactions
	dXdt[c12] = A[c12] * ( Y[he4] * Y[he4] * Y[he4] * rate[r3a]
				- Y[c12] * rate[rg3a]
				+ Y[o16] * rate[roga]
				- Y[c12] * Y[he4] * rate[rcag]
				- 2.0 * Y[c12] * Y[c12] * rate[r1212]
				- Y[c12] * Y[o16] * rate[r1216]);
    

	// 16O reactions
	dXdt[o16] = A[o16] * (- Y[o16] * rate[roga]
				+ Y[c12] * Y[he4] * rate[rcag]
				- Y[c12] * Y[o16] * rate[r1216]
				- 2.0 * Y[o16] * Y[o16] * rate[r1616]
				- Y[o16] * Y[he4] * rate[roag]
				+ Y[ne20] * rate[rnega] );

	// 20Ne reactions
	dXdt[ne20] = A[ne20] * (Y[c12] * Y[c12] * rate[r1212]
				+ Y[o16] * Y[he4] * rate[roag]
				- Y[ne20] * rate[rnega]
				+ Y[mg24] * rate[rmgga]
				- Y[ne20] * Y[he4] * rate[rneag]);

	// 24Mg reactions
	dXdt[mg24] = A[mg24] * (0.5 * Y[c12] * Y[o16] * rate[r1216]
				- Y[mg24] * rate[rmgga]
				+ Y[ne20] * Y[he4] * rate[rneag]
				+ Y[si28] * rate[rsiga]
				- Y[mg24] * Y[he4] * rate[rmgag] ) ;

	// 28Si reactions
	dXdt[si28] = A[si28] * (0.5 * Y[c12] * Y[o16] * rate[r1216]
                + Y[o16] * Y[o16] * rate[r1616]
				- Y[si28] * rate[rsiga]
				+ Y[mg24] * Y[he4] * rate[rmgag]
				- Y[he4] * rate[rsi2ni]
				+ Y[ni56] * rate[rni2si]);

	// 56Ni reactions
	dXdt[ni56] = A[ni56] * ( Y[he4] * rate[rsi2ni]
				- Y[ni56] * rate[rni2si] );

    //Cleaning the house
    Y.clear();
}


void Nuclear_Network::Calculate_Jacobian(const double &rho,
 const double &T, const std::vector<double> &X, std::vector< std::vector<double> > &dfdy)
{
    std::vector<double> Y(7);
    const double A[7] = {4.0,12.0,16.0,20.0,24.0,28.0,56.0};
    
    for (unsigned int i = 0; i < 7; i++)
        Y[i] = X[i]/A[i];
    
    Calculate_Rates(rho,T,Y);
    
  //#pragma omp parallel for
    for (unsigned int i = 0; i < 7; i++) {
        for (unsigned int j = 0; j < 7; j++)
            dfdy[i][j] = 0.0;
    }
    
    // 4He jacobian elements
    dfdy[he4][he4] = -9.0 * Y[he4] * Y[he4] * rate[r3a]
				- Y[c12] * rate[rcag]
				- Y[o16] * rate[roag]
				- Y[ne20] * rate[rneag]
				- Y[mg24] * rate[rmgag]
				- 7.0 * rate[rsi2ni]
				- 7.0 * rsi2nida * Y[he4]
                + 7.0 * rni2sida * Y[ni56];
    
    dfdy[he4][c12] = (A[he4]/A[c12])*(3.0 * rate[rg3a]
				- Y[he4] * rate[rcag]
				+ 2.0 * Y[c12] * rate[r1212]
				+ 0.5 * Y[o16] * rate[r1216]);
    
    dfdy[he4][o16] = (A[he4]/A[o16])*(rate[roga]
				+ 0.5 * Y[c12] * rate[r1216]
				+ 2.0 * Y[o16] * rate[r1616]
				- Y[he4] * rate[roag]);
    
    dfdy[he4][ne20] = (A[he4]/A[ne20]) * (rate[rnega] - Y[he4] * rate[rneag]);
    
    dfdy[he4][mg24] = (A[he4]/A[mg24]) * ( rate[rmgga] - Y[he4] * rate[rmgag]);
    
    dfdy[he4][si28] = (A[he4]/A[si28]) * ( rate[rsiga]
				- 7.0 * rsi2nidsi * Y[he4]);
    
    dfdy[he4][ni56] = (A[he4]/A[ni56]) * (7.0 * rni2si);
    
    // 12C jacobian elements
    dfdy[c12][he4] = (A[c12]/A[he4]) * (3.0 * Y[he4] * Y[he4] * rate[r3a]
				- Y[c12] * rate[rcag]);
    
    dfdy[c12][c12] = - rate[rg3a]
				- Y[he4] * rate[rcag]
				- 4.0 * Y[c12] * rate[r1212]
				- Y[o16] * rate[r1216];
    
    dfdy[c12][o16] = (A[c12]/A[o16]) * ( rate[roga]
				- Y[c12] * rate[r1216]);
    
    
    // 16O reactions
    dfdy[o16][he4] = (A[o16]/A[he4]) * (Y[c12] * rate[rcag]
				- Y[o16] * rate[roag]);
    
    dfdy[o16][c12] = (A[o16]/A[c12]) * (Y[he4] * rate[rcag]
				- Y[o16] * rate[r1216]);
    
    dfdy[o16][o16] = (A[o16]/A[o16]) * (- rate[roga]
				- Y[c12] * rate[r1216]
				- 4.0 * Y[o16] * rate[r1616]
				- Y[he4] * rate[roag]);
    
    dfdy[o16][ne20] = (A[o16]/A[ne20]) * (rate[rnega]);
    
    // 20Ne jacbian elements
    dfdy[ne20][he4] = (A[ne20]/A[he4]) * (Y[o16] * rate[roag]
				- Y[ne20] * rate[rneag]);
    
    dfdy[ne20][c12] = (A[ne20]/A[c12]) * (2.0 * Y[c12] * rate[r1212]);
    
    dfdy[ne20][o16] = (A[ne20]/A[o16]) * (Y[he4] * rate[roag]);
    
    dfdy[ne20][ne20] = - rate[rnega] - Y[he4] * rate[rneag];
    
    dfdy[ne20][mg24] = (A[ne20]/A[mg24]) * (-rate[rmgga]);
    
    // 24Mg jacobian elements
    dfdy[mg24][he4] = (A[mg24]/A[he4]) * (Y[ne20] * rate[rneag]
				- Y[mg24] * rate[rmgag]);
    
    dfdy[mg24][c12] = (A[mg24]/A[c12]) * (0.5 * Y[o16] * rate[r1216]);
    
    dfdy[mg24][o16] = (A[mg24]/A[o16]) * (0.5 * Y[c12]  * rate[r1216]);
    
    dfdy[mg24][ne20] = (A[mg24]/A[ne20]) * (Y[he4] * rate[rneag]);
    
    dfdy[mg24][mg24] = - rate[rmgga]
				- Y[he4] * rate[rmgag];
    
    dfdy[mg24][si28] = (A[mg24]/A[si28]) * rate[rsiga];
    
    // 28Si jacobian elements
    dfdy[si28][he4] = (A[si28]/A[he4]) * (Y[mg24] * rate[rmgag]
				- rate[rsi2ni]
				- rsi2nida * Y[he4] + rni2sida * Y[ni56]);
    
    dfdy[si28][c12] = (A[si28]/A[c12]) * (0.5 * Y[o16] * rate[r1216]);
    
    dfdy[si28][o16] = (A[si28]/A[o16]) * (0.5 * Y[c12] * rate[r1216]
                + 2.0 * Y[o16] * rate[r1616]);
    
    dfdy[si28][mg24] = (A[si28]/A[mg24]) * (Y[he4] * rate[rmgag]);
    
    dfdy[si28][si28] = - rate[rsiga]
				- rsi2nidsi * Y[he4];
    
    dfdy[si28][ni56] = (A[si28]/A[ni56]) * rate[rni2si];
    
    // 56Ni jacobian elements
    dfdy[ni56][he4] =  (A[ni56]/A[he4]) * (rate[rsi2ni] + rsi2nida * Y[he4]
                - rni2sida * Y[ni56]);
    
    dfdy[ni56][si28] = (A[ni56]/A[si28]) * rsi2nidsi * Y[he4];
    
    dfdy[ni56][ni56] = -rate[rni2si];
    
    //Cleaning the house
    Y.clear();
}

// All the reaction rates below are from Caughlan and Fowler 1988, except ca40(ag)
// and its inverse, which is from Woosley et al. 1978. One may wish to substitute
// the corresponding reaction rates from the NACRE (Angulo et al. 1999) or
// Rauscher & Thielemann (1998) libraries

void Nuclear_Network::Calculate_Rates(const double &rhop, const double &Tp,
                                      const std::vector<double> &Y){
    
    rho = rhop;
    T = Tp;

    rate.resize(17);
    
    // Some temperature factors
    // limit T9 to 10, except for the reverse ratios
    
    TT9 = T*1.e-9;
    T9r = TT9;
    T9 = MIN(TT9,10.0);
    T912 = sqrt(T9);
    T913 = pow(T9,1/3.);
    T923 = Q(T913);
    T943 = T9*T913;
    T953 = T9*T923;
    T932 = T9*T912;
    T92 = Q(T9);
    T93 = T9*T92;
    T972 = T912*T93;
    T9r32 = T9r*sqrt(T9r);
    
    T9i = 1./T9;
    T9i13 = 1./T913;
    T9i23 = 1./T923;
    T9i32 = 1./T932;
    T9i12 = 1./T912;
    T9ri = 1./T9r;
    
    TripleAlpha_Rate();
    CarbonFusion_Rate();
    CarbonOxygen_Rate();
    OxygenFusion_Rate();
    CarbonAlpha_Rate();
    OxygenAlpha_Rate();
    NeonAlpha_Rate(); 
    MagnesiumAlpha_Rate(); 
    CalciumAlpha_Rate(); 
    Silicon2Nickel_Rate(Y);

    //Electron_Screening(Y);
}

struct root_He {
    double ratio1, ratio2;

    double operator()(double y){

        return (pow(y,14.0)/ratio1 + y + pow(y,7.0)/ratio2- 1.0);
    }
};

void Nuclear_Network::NuclearStatisticalEquilibrium(const double &rhop, const double &Tp, std::vector<double> &X){
    
    using PhysConstants::k_B;
    using PhysConstants::N_A;
    using PhysConstants::MeV2erg;
    
    //Binding Energy in keV
    const double B[7] = {28295.673,92161.753,127619.336,160644.852,198256.887,236536.884,483987.827};
    
    for (unsigned int i = 0; i < X.size(); i++) {
        X[i] = 0.0; //All abundances are zero, except He, Si, and Ni
    }
    
    const double theta = 5.943e33*pow(Tp/1.0e9,1.5);
    
    //Abundance of He over Ni, following the NSE (Saha equation)
    double HeoverNi = 5.03071e16*pow(theta/(rhop*N_A),13.0)*exp((14*B[0]-B[6])*MeV2erg*1.e-3/(k_B*Tp));
    
    //Abundance of Si over Ni, following the NSE (Saha equation)
    double HeoverSi = 8.28237e6*pow(theta/(rhop*N_A),6.0)*exp((7*B[0]-B[5])*MeV2erg*1.e-3/(k_B*Tp));
    
    //The sum of abundances must be one!
    root_He func;
    
    func.ratio1 = HeoverNi;
    func.ratio2 = HeoverSi;
    
    X[0] = zbrent(func,0.0,1.0,1e-16);
    X[5] = pow(X[0],7.0)/HeoverSi;
    X[6] = (1.0 - X[5] - X[0]);
}

double ScreeningFactor(const double &rhop, const double &T,
                       const double &Z1, const double &Z2,const double &Ye){
	
	double T6 = T*1.e-6;

	return (1+0.205*(pow(Z1+Z2,5/3.)-pow(Z1,5/3.)-pow(Z2,5/3.))*pow(Ye*rhop,1/3.)/T6);
}

void Nuclear_Network::Electron_Screening(const std::vector<double> &Y){

    double Z[7] = {2,6,8,10,12,14,28};

    double Ye = Z[0]*Y[0] + Z[1]*Y[1] + Z[2]*Y[2]+ Z[3]*Y[3] + Z[4]*Y[4] + Z[5]*Y[5] + Z[6]*Y[6]; 

    rate[r3a]*= ScreeningFactor(rho,T,Z[0],Z[0],Ye);

    rate[r1212]*= ScreeningFactor(rho,T,Z[1],Z[1],Ye);

    rate[r1216] *= ScreeningFactor(rho,T,Z[1],Z[2],Ye);

    rate[r1616] *= ScreeningFactor(rho,T,Z[2],Z[2],Ye);

    rate[rcag] *= ScreeningFactor(rho,T,Z[1],Z[0],Ye);

    rate[roag] *= ScreeningFactor(rho,T,Z[2],Z[0],Ye);

    rate[rneag] *= ScreeningFactor(rho,T,Z[3],Z[0],Ye);

    rate[rmgag] *= ScreeningFactor(rho,T,Z[4],Z[0],Ye);

    rate[rcaag] *= ScreeningFactor(rho,T,20,Z[0],Ye);
}


// triple alpha to 12C reaction
void Nuclear_Network::TripleAlpha_Rate()
{
	double r2abe = 7.40e5 * T9i32 * exp(-1.0663*T9i) +
					4.164e9 * T9i23 * exp(-13.49*T9i13 - T92/0.009604) *
					(1. + 0.031*T913 + 8.009*T923 + 1.732*T9 + 
					49.883*T943 + 27.426*T953);

	double rbeac = 130. * T9i32 * exp(-3.3364*T9i) +
					2.510e7 * T9i23 * exp(-23.57*T9i13 - T92/0.055225) * 
					(1. + 0.018*T913 + 5.249*T923 + 0.650*T9 +
					19.176*T943 + 6.034*T953);

	double term;		

	if (T9 > 0.08) term = 2.90e-16 * (r2abe*rbeac) +
						  0.1 * 1.35e-7 * T9i32 * exp(-24.811*T9i);

	else term = 2.90e-16 * (r2abe*rbeac) * 
				(0.01 + 0.2*(1. + 4.*exp(-pow(0.025*T9i,3.263)))/
				(1. + 4.*exp(-pow(T9/0.025,9.227)))) + 
				0.1 * 1.35e-7 * T9i32 *exp (-24.811*T9i);

	rate[r3a] = term * (rho*rho)/6.;

	double rev = 2.0e20 * T9r32 * exp(-84.424*T9ri);

	rate[rg3a] = rev * term;
}

// 12C + 12C reaction
void Nuclear_Network::CarbonFusion_Rate()
{
	double T9a = T9/(1. + 0.0396*T9);
	double T9a13 = pow(T9a,1/3.);
	double T9a56 = pow(T9a,5/6.);

	double term = 4.27e26 * T9a56 * T9i32 * 
					exp(-84.165/T9a13 - 2.12e-3*T9*T9*T9);

	rate[r1212] = 0.5 * rho * term;

}

// 12C + 16O reaction
void Nuclear_Network::CarbonOxygen_Rate()
{
	if (T9 >= 0.5){
		double T9a = T9/(1. + 0.055*T9);
		double T9a13 = pow(T9a,1/3.);
		double T9a23 = T9a13*T9a13;
		double T9a56 = pow(T9a,5/6.);

		double term = 1.72e31 * T9a56 * T9i32 * exp(-106.594/T9a13) /
						(exp(-0.18*T9a*T9a) + 1.06e-3*exp(2.562*T9a23));

		rate[r1216] = rho * term;
        
	}

	else rate[r1216] = 0.;
}

// 16O + 16O reaction
void Nuclear_Network::OxygenFusion_Rate()
{
	double term = 7.10e36 * T9i23 * 
					exp(-135.93*T9i13 - 0.629*T923 -
						0.445*T943 + 0.0103*T92);

	rate[r1616] = 0.5 * rho * term;
}

// 12C(a,g)16O and inverse
void Nuclear_Network::CarbonAlpha_Rate()
{
	double term = 1.04e8 * exp(-32.120*T9i13 - T92/12.222016) / (T92 * Q(1. + 0.0489*T9i23))
                    + 1.76e8 * exp(-32.120*T9i13)/(T92*Q(1. + 0.2654*T9i23))
					+ 1.25e3*T9i32*exp(-27.499*T9i)
					+ 1.43e-2*T92*T93*exp(-15.541*T9i);

	term = 1.7 * term;

	rate[rcag] = rho * term;

	double rev = 5.13e10 * T9r32 * exp(-83.111*T9ri);

	rate[roga] = rev * term;
}

// 16O(a,g)20Ne and inverse
void Nuclear_Network::OxygenAlpha_Rate()
{
	double term1 = 9.37e9 * T9i23 *
					exp(-39.757*T9i13 - T92/2.515396);

	double term2 = 62.1 * T9i32 * exp(-10.297*T9i) +
					538.*T9i32*exp(-12.226*T9i) +
					13.*T92*exp(-20.093*T9i);

	double term = term1 + term2;

	rate[roag] = rho * term;

	double rev = 5.65e10 * T9r32 * exp(-54.937*T9ri);

	rate[rnega] = rev * term;
}

// 20Ne(a,g)24Mg and inverse
void Nuclear_Network::NeonAlpha_Rate()
{
	double term1 = 4.11e11* T9i23 * 
					exp(-46.766*T9i13 - T92/4.923961) *
					(1.0 + 0.009* T913 + 0.882*T923 + 0.055*T9 +
					0.749*T943 + 0.119*T953);

	double term2 = 5.27e3 * T9i32 * exp(-15.869*T9i) + 
					6.51e3 * T912 * exp(-16.223*T9i);

	double term3 = 0.1 * (42.1 * T9i32 * exp(-9.115*T9i) +
					32.0 * T9i23 * exp(-9.383*T9i));

	double term = (term1+term2+term3)/(1. + 5.*exp(-18.960*T9i));

	rate[rneag] = rho * term;

	double rev = 6.01e10 * T9r32 * exp(-108.059*T9ri);

	rate[rmgga] = rev * term;
}

// 24Mg(a,g)28Si and inverse
void Nuclear_Network::MagnesiumAlpha_Rate()
{
	double term1 = 47.8 * T9i32 * exp(-13.506*T9i) +
					2.38e3 * T9i32 * exp(-15.218*T9i) +
					2.47e2 * T932 * exp(-15.147*T9i) + 
					0.1 * (1.72e-9 * T9i32 * exp(-5.028*T9i) +
						1.25e-3 * T9i32 * exp(-7.929*T9i) + 
						24.3 * T9i *exp(-11.523*T9i));

	double term = term1/(1. + 5.*exp(-15.882*T9i));

	rate[rmgag] = rho * term;

	double rev = 6.27e10 * T9r32 * exp(-115.862*T9ri);

	rate[rsiga] = rev * term;
}

// 40Ca(a,g)44Ti and inverse
void Nuclear_Network::CalciumAlpha_Rate()
{
	double term = 4.66e24 * T9i23 * exp(-76.435 * T9i13 *
					(1. + 1.650e-2 * T9 + 5.973e-3*T9*T9 -
						3.889e-4*T9*T9*T9));

	rate[rcaag] = rho * term;

	double rev = 6.843e10 * T9r32 * exp(-59.510*T9ri);

	rate[rtiga] = rev * term;
}

// 28Si to 56Ni and inverse
void Nuclear_Network::Silicon2Nickel_Rate(const std::vector<double> &Y)
{
	if (T9 > 2.5 && Y[c12] + Y[o16] <= 0.004){
		rate[rsi2ni] = pow(T9i32,3) * exp(239.42*T9i-74.741) *
					pow(rho,3) * pow(Y[he4],3) * rate[rcaag]*Y[si28];

		rate[rni2si] = MIN(1.e20,pow(T932,3) * exp(-274.12*T9i+74.914) *
					 rate[rtiga]/(pow(rho,3)* pow(Y[he4],3)) );
        
        rsi2nida = 3. * pow(T9i32,3) * exp(239.42*T9i-74.741) *
                pow(rho,3) * pow(Y[he4],2) * rate[rcaag]*Y[si28];
        rsi2nidsi = pow(T9i32,3) * exp(239.42*T9i-74.741) *
                pow(rho,3) * pow(Y[he4],3) * rate[rcaag];

        if (rate[rni2si] == 1.e20 ) rni2sida = 0.;
        else rni2sida = -3. * rate[rni2si]/ Y[he4] ;

	}

	else {
		rate[rsi2ni] = 0.;
		rate[rni2si] = 0.;
        rsi2nida = 0.;
        rsi2nidsi = 0.;
        rni2sida = 0.;
	}
}

double Nuclear_Network::EnergyGenerated(std::vector<double> &dXdt){
    
    using PhysConstants::N_A;
    using PhysConstants::MeV2erg;
    
    const double B[7] = {28295.673,92161.753,127619.336,160644.852,198256.887,236536.884,483987.827};
    const double A[7] = {4.0,12.0,16.0,20.0,24.0,28.0,56.0};
    
    double Eps = 0.;
    for (unsigned int i = 0; i < dXdt.size(); i++) {
        Eps += N_A*B[i]*MeV2erg*1.e-3*dXdt[i]/A[i];
    }
    
    return Eps;
}

double Nuclear_Network::EnergyGeneratedIntegrated(const std::vector<double> &X,
		 const std::vector<double> &Xi, const double &dt){
    
    using PhysConstants::N_A;
    using PhysConstants::MeV2erg;
    
    const double B[7] = {28295.673,92161.753,127619.336,160644.852,198256.887,236536.884,483987.827};
    const double A[7] = {4.0,12.0,16.0,20.0,24.0,28.0,56.0};
    
    double Eps = 0.;
    for (unsigned int i = 0; i < X.size(); i++) {
        Eps += -N_A*B[i]*MeV2erg*1.e-3*(Xi[i]-X[i])/(A[i]*dt);
    }
    
    return Eps;
}
