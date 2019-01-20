#ifndef _PHOTONS_HPP_
#define _PHOTONS_HPP_

#include <iostream>
#include <vector>
#include <stdlib.h>

#include "../cpp_constants_CGS.hpp"

struct Photons {
public:
    double n, p, u, s, dpdrho, dpdT, dudrho, dudT, dsdrho, dsdT;
    void Evaluate(const double& rho_in, const double& T_in);

private:
    double rho, T;
    double EnergyDensity();

    double Pressure();

    double EntropyDensity();
};

//=========================================
//  PHOTON GAS
// =========================================

void Photons::Evaluate(const double& rho_in, const double& T_in)
{
    rho = rho_in;
    T = T_in;

    p = Pressure();
    u = EnergyDensity()/rho; //specific internal energy
    s = EntropyDensity()/rho; //specific entropy
    dpdrho = 0.0 ;
    dpdT = 4.*p/T;
    dudrho = -u/rho ; 
    dudT = 4.*u/T;
    dsdrho = - s /rho;
    dsdT = 3.*s/T;
}

double Photons::Pressure()
{
    using PhysConstants::a;
 
    return a*pow(T,4.)/3.;
}

double Photons::EnergyDensity()
{       
    return 3.*p;
}

double Photons::EntropyDensity()
{   
    return (rho * u + p) / T;
}

#endif