#ifndef _IONS_HPP_
#define _IONS_HPP_

#include <iostream>
#include <vector>
#include <stdlib.h>

#include "../cpp_constants_CGS.hpp"

struct Ions {
public:
    double n, p, u, s, dpdrho, dpdT, dudrho, dudT, dsdrho, dsdT;
    void Evaluate(const double& rho_in, const double& T_in, const double &Yi_in);

private:
    double rho, T, Yi;
    double EnergyDensity();

    double Pressure();

    double EntropyDensity();

};

//=========================================
//  ION GAS
// =========================================

void Ions::Evaluate(const double& rho_in, const double& T_in, const double &Yi_in){

    using PhysConstants::N_A;
    using PhysConstants::k_B;

    rho = rho_in;
    T = T_in;
    Yi = Yi_in;
    
    n = rho*N_A*Yi; // number of free electrons

    p = Pressure();
    u = EnergyDensity()/rho; //specific internal energy
    s = EntropyDensity()/rho; //specific entropy
    dpdrho = N_A*Yi * k_B * T ;
    dpdT = n * k_B;
    dudrho = 0. ; 
    dudT = 1.5* n * k_B/rho;
    dsdrho = - k_B * N_A * Yi /rho;
    dsdT = 1.5 * k_B * N_A *Yi/T;
}

double Ions::Pressure()
{
    using PhysConstants::k_B;
 
    return n * k_B * T;
}

double Ions::EnergyDensity()
{       
    return 1.5 * p;
}

double Ions::EntropyDensity()
{
    using PhysConstants::k_B;
    using PhysConstants::pi;
    using PhysConstants::h;
    using PhysConstants::amu;

    double m = amu/Yi;
    
    return n*k_B*(log(pow(4*pi*m*u/(3*n*Q(h)),1.5)/n) + 2.5);
}

#endif