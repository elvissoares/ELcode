/*Stellar_Structure_Equations

    General Description:
    ====================
 
        This header file contains the function interface definitions for the basic equations of 
        stellar structure.  The file also contains a driver function interface definition that selects 
        among the required equations for the Runge Kutta routine.
---------------------------------------------------------------------*/

#ifndef _STRUCTUREEQNS_HPP_
#define _STRUCTUREEQNS_HPP_

#include <iostream>
#include <cmath>
#include <vector>

#include "cpp_constants_CGS.hpp"

//=====================================================================
// Parameter \[Ksi] = Subscript[R, i - 1]/Subscript[R, i]
//=====================================================================
double Ksi(const double &Rout, const double &Rin)
{
    return Rin / Rout;
};

//=====================================================================
// Função que calcula a energia gravitacional Vg total da estrela
// para uma dada configuração {R_i}
//=====================================================================

double f(const double& ksi)
{
    double denom = Q(ksi)+ksi+1.;
	return (2. + 4. * ksi + 6. * Q(ksi) + 3. * C(ksi)) / Q(denom);
};

double g(const double& ksi)
{
    double denom = Q(ksi)+ksi+1.;
	return (5. * (1. + ksi) / denom);
};


//=====================================================================
// Time derivate of parameter Ksi(R,j) \[Ksi(R,j)]'= d\[Ksi(R,j)]/dt
//=====================================================================

double Ksidot(const double &Rout, const double &Rin, const double &Rdot_out, const double &Rdot_in)
{
    return (Rout * Rdot_in - Rin * Rdot_out) / Q(Rout);
};

//=====================================================================
// Matriz de energia cinética (T)
//=====================================================================

double T_11(const double& ksi, const double& m)
{
	const double threefifths = 0.6;
	double denom = 1.+ksi+Q(ksi); //(1 + \[Ksi] + \[Ksi]^2)^3
	
	return (threefifths*C(ksi)*(5. +6.*ksi +3.*Q(ksi) +C(ksi)) * m / C(denom));
};

double T_12(const double& ksi, const double& m)
{
	const double ninetenths = 0.9;
	double denom = 1.+ksi+Q(ksi);
	
	return (ninetenths*Q(ksi)*(1. +3.*ksi +Q(ksi)) * m / C(denom));
};

double T_22(const double& ksi, const double& m)
{
	const double threefifths = 0.6;
	double denom = 1.+ksi+Q(ksi);
	
	return (threefifths*(1. +3.*ksi +6.*Q(ksi) +5.*C(ksi)) * m / C(denom));
};

//=====================================================================
// Derivada temporal da Matriz de energia cinética (Q=-(1/2)*dT/dt)
//=====================================================================

double Q_11(const double& ksi, const double& ksidot, const double& m)
{
	const double ninetenths = 0.9;
	double denom = 1.+ksi+Q(ksi); // 1 + \[Ksi] + \[Ksi]^2
	
	return (-ninetenths*Q(ksi)*(5. +8.*ksi +2.*Q(ksi))*ksidot*m/pow(denom,4));
};

double Q_12(const double& ksi, const double& ksidot, const double& m)
{
	const double ninetenths = 0.9;
	double denom = 1.+ksi+Q(ksi);
	
	return (ninetenths*ksi*(-1. -4.*ksi +4.*C(ksi) +pow(ksi,4.))*ksidot*m/pow(denom,4));
};

double Q_22(const double& ksi, const double& ksidot, const double& m)
{
	const double ninetenths = 0.9;
	double denom = 1.+ksi+Q(ksi);
	
	return (ninetenths*Q(ksi)*(2. +8.*ksi +5.*Q(ksi))*ksidot*m/pow(denom,4));
};

//=====================================================================
// funções para o Termo de Força (f, g)
//=====================================================================

double f_1(const double& ksi)
{
	const double threefifths = 0.6;
	double denom = 1.+ksi+Q(ksi); // 1 + \[Ksi] + \[Ksi]^2
	
	return (threefifths*(1. +3.*ksi +6.*Q(ksi) +5.*C(ksi))/C(denom));
};

double f_2(const double& ksi)
{
	const double ninetenths = 0.9;
	double denom = 1.+ksi+Q(ksi);
	
	return (ninetenths*pow(ksi,4)*(1. +3.*ksi +Q(ksi))/C(denom));
};

double g_1(const double& ksi)
{
	const double threehalfs = 1.5;
	double denom = 1.+ksi+Q(ksi); // 1 + \[Ksi] + \[Ksi]^2
	
	return (threehalfs*(1. +2.*ksi)/Q(denom));
};

double g_2(const double& ksi)
{
	const double threehalfs = 1.5;
	double denom = 1.+ksi+Q(ksi);
	
	return (threehalfs*C(ksi)*(2. +ksi)/Q(denom));
};

#endif
