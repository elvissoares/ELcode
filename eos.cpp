#include "eos.hpp"
#include "cpp_constants_CGS.hpp"

#include <cmath>
#include <stdio.h>
#include <string>
#include <fstream>

#include "util/roots.hpp"

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

//=========================================
//  EOS Struct
// =========================================

void EoS::Evaluate(const double& rho, const double& T, const std::vector<double> &Y)
{

    // 0 for pressuless EoS
    if (id_eos == 0) {
        P = 0.;
        dPdrho = 0.;
        dPdT = 0.;
        
        u = 0.;
        dudrho = 0.;
        dudT = 0.;
        
        s = 0.;
        dsdrho = 0.;
        dsdT = 0.;
        
        gamma = 0.;
        cs = 0.;
        cv = 0.;
        cp = 0.;
    }
    
    // 1 for polytropic EoS
    if (id_eos == 1) {

        using PhysConstants::N_A;
        using PhysConstants::k_B;

        const double R = N_A*k_B;
        
        gamma = 7/5.;
        double K = 1.0 + (R/pow(rho,gamma-1))*(T-1e7);
        
        P = K*pow(rho,gamma); 

        cv = R/(gamma-1.);
        cp = gamma*cv;

        u = (P/rho)/(gamma-1);
        dudrho = P/Q(rho);
    }
    
    // 2 for ions + photons + electrons Helmholtz EoS
    if (id_eos == 2){
        
        Yi = Ion_MolarAbundance(Y);
        Zi = Ion_MeanCharge(Y);
        Ye = ElectronsPerBaryons(Y);
        
        ions.Evaluate(rho,T,Yi);
        photons.Evaluate(rho,T);
        pairs.Evaluate(rho,T,Ye);
        
        P = ions.p + photons.p + pairs.p;
        dPdrho = ions.dpdrho + photons.dpdrho + pairs.dpdrho;
        dPdT = ions.dpdT + photons.dpdT + pairs.dpdT;
        
        u = ions.u + photons.u + pairs.u;
        dudrho = ions.dudrho + photons.dudrho + pairs.dudrho;
        dudT = ions.dudT + photons.dudT + pairs.dudT;
        
        s = ions.s + photons.s + pairs.s;
        dsdrho = ions.dsdrho + photons.dsdrho + pairs.dsdrho;
        dsdT = ions.dsdT + photons.dsdT + pairs.dsdT;
        
        gamma = (rho/P)*dPdrho;
        cs = sqrt((dPdrho - dPdT*dsdrho/dsdT));
        cv = dudT;
        cp = dudT + dPdT/rho;
    }
    
}

//Calculate the mean molecular weight of the gas
double EoS::Ion_MolarAbundance(const std::vector<double> &X)
{
    double A[7] = {4,12,16,20,24,28,56};

    double sum = 0.;
    for (std::vector<int>::size_type i = 0; i < X.size(); i++)
        sum += X[i]/A[i];
    
    return sum;
}

//Calculate the mean molecular weight of the gas
double EoS::Ion_MeanCharge(const std::vector<double> &X)
{
    double Z[7] = {2,6,8,10,12,14,28};
    
    double sum = 0.;
    for (std::vector<int>::size_type i = 0; i < X.size(); i++)
        sum += Z[i]*X[i];
    
    return sum;
}

double EoS::ElectronsPerBaryons(const std::vector<double> &X)
{
    //For complete ionization
    const double Z[7] = {2,6,8,10,12,14,28};
    double A[7] = {4,12,16,20,24,28,56};
    
    double sum = 0.;
    for (std::vector<int>::size_type i = 0; i < X.size(); i++)
        sum += Z[i]*X[i]/A[i];
    
    return sum;
}