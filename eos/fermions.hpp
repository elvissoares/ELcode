#include "fermi.hpp"

#ifndef EOSPAIRS_HPP_
#define EOSPAIRS_HPP_

#include <cmath>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

struct Fermions: public FermiDirac, FermiDiracDBeta, FermiDiracDEta{
 private:
    double beta, T, rho, Ye;
    double F32, F52, dF32deta, dF52deta, dF32dbeta, dF52dbeta;
    double F12_n, F32_n, dF12deta_n, dF32deta_n, dF12dbeta_n, dF32dbeta_n;
    
    void EvaluateIntegrals(const double& eta, const double& beta);
    double Pressure();
    double EnergyDensity();
    double dPdRho();
    double dPdT();
    double dUdRho();
    double dUdT();

    double Free(const double &rhoin, const double &Tin);
    double d2PdRhodT();
    double d2PdT2();
    double d3PdT2dRho();

    //friend struct root_Eta;
    double Eta_Search(const double& n, const double& T);

 public: 
    double B_mu, A;
 	double eta, n, p, u, s, dpdrho, dpdT, dudrho, dudT, dsdrho, dsdT;
    double d3fdT2drho,d3fdrho2dT,d4fdrho2dT2 ;
    Fermions();
    void Evaluate(const double& rho, const double& T, const double &Ye);
    void Evaluate_FreeEnergyDerivatives();
    double Beta(const double& );
    
};

struct InterpolatedFermions {
public:
    double p, u, s, dpdrho, dpdT, dudrho, dudT, dsdrho, dsdT;
    double Rhomin, Rhomax, Tmin, Tmax;
    void Evaluate(const double &Rhoin, const double &Tin, const double &Ye);

    InterpolatedFermions(){
        Read_Table();
    }

private:
    std::vector<double> Rho, T;
    std::vector< std::vector<double> > f, ft, fd, ftt, fdd, fdt, fddt, fdtt, fddtt;
    double psi0(const double &z);
    double dpsi0(const double &z);
    double ddpsi0(const double &z);

    double psi1(const double &z);
    double dpsi1(const double &z);
    double ddpsi1(const double &z);

    double psi2(const double &z);
    double dpsi2(const double &z);
    double ddpsi2(const double &z);

    double herm5(const int &i, const int &j, const double &w0t,
    const double &w1t, const double &w2t, const double &w0mt, const double &w1mt, 
    const double &w2mt, const double &w0d, const double &w1d, const double &w2d,
    const double &w0md, const double &w1md, const double &w2md);

    void Read_Table();
};

#endif