//=====================================================================
// Author: Elvis Soares
// Date: June 2015
// last modification: May 2016
//=====================================================================

#ifndef _EOS_HPP_
#define _EOS_HPP_

#include <iostream>
#include <vector>
#include <stdlib.h>

#include "eos/fermions.hpp"

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

struct EoS {
    
public:
    InterpolatedFermions pairs;
    Ions ions;
    Photons photons;

    double Yi, Zi, Ye, P, u, s, dPdrho, dPdT, dudrho, dudT, dsdrho, dsdT, gamma, cs, cv, cp;
    
    unsigned int id_eos;
    
    EoS(unsigned int ideos = 2): id_eos(ideos) {}
    
    void Evaluate(const double& rho_in, const double& T_in, const std::vector<double> &Y_in);

    double Ion_MolarAbundance(const std::vector<double> &Y);

    double Ion_MeanCharge(const std::vector<double> &Y);
    
    double ElectronsPerBaryons(const std::vector<double> &X);

};

#endif