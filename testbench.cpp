#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cassert>
#include <vector>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 0
#define omp_get_thread_num() 0
#endif

#include "cpp_constants_CGS.hpp"
#include "star.hpp"
#include "user_IO.hpp"
#include "util/nr.hpp"
#include "util/mins_ndim.hpp"
#include "util/odeint.hpp"
#include "eos/fermi.hpp"
#include "thermonuclear.hpp"
#include "util/roots.hpp"
#include "conductivity.hpp"

#define tab "\t"

#define test1 0 // Produz uma tabela da energia livre de Helmholtz em funçao de log rho e log T
#define test1b 0 // Phase diagram
#define test2 1// Calcula as qtds termodinâmicas para diferentes temperaturas
#define test3 0 // Teste para Nuclear Network iso7
#define test4 0 // Taxa de reação para He4
#define test5 0 // Teste para condutividade térmica por elétrons
#define test10 0 //Test to evolve white dwarf during 1e8 years
#define test11 0 //Test to evolve many white dwarfs during 1e8 years
#define test12 0 //Free-fall test (Ok!)
#define test13 0 //Sod Shock Wave test (Ok!)
#define test14 0 //Sedov's Blast Wave test
#define test20 0
#define test21 0//Dynamics with artificial mass pertubation

using namespace std;

template <typename T>
std::string to_string(T const& value) {
    stringstream sstr;
    sstr << value;
    return sstr.str();
}

struct rootHe: public Nuclear_Network {
    double Rho, tau; // n-point and T-point
    std::vector<double> Y;
    EoS eos;

    double operator()(double logT){

        double T = pow(10,logT);

        std::vector<double> dYdt(7);

        Calculate_Derivatives(Rho,T,Y,dYdt);

        double eps = EnergyGenerated(dYdt);

        eos.Evaluate(Rho,T,Y);

        double taup = fabs(eos.cp) * T / eps;

        dYdt.clear();
    
        return (tau-taup);
    }
};

struct Network: public Nuclear_Network {
    double Rho, T;
    
    void operator() (const double &t, const std::vector<double> &y, std::vector<double> &dydt)
    {
        Calculate_Derivatives(Rho,T,y,dydt);
    }
    
    void jacobian(const double &t, const std::vector<double> &y, std::vector<double> &dydt, std::vector< std::vector<double> > &dfdy)
    {
        Calculate_Jacobian(Rho,T,y,dfdy);
    }
};

int main (void) {
    
    std::cerr << "The max number of threads is " << omp_get_max_threads() << std::endl;
   
   	time_t begin, end;

    begin = time(0);

#if test1 // Produz uma tabela da energia livre de Helmholtz em funçao de log rho e log T
{
    double logrho_lo = -6, logrho_hi = 11;
    double drho = 0.1;
    
    double logT_lo = 3, logT_hi = 11;
    double dT = 0.1;
    
    string OutName = "output/eos/free_energy-pairs.dat";
    ofstream OutFile(OutName.c_str());
    
    OutFile << "# EoS grid of quantities for by Rho/mu_e and T "<< endl;
    OutFile << "# Number of x1 points =" << (logrho_hi-logrho_lo)/drho +1 << endl;
    OutFile << "# Number of x2 points =" << (logT_hi-logT_lo)/dT +1 << endl;
    OutFile << "# Rho/mu_e  T  F  dFdrho  dFdT  d2FdRho2  d2FdT2   d2FdRhodT" << std::endl;
    
    Fermions pairs;
    
    double Rho, T, F, dFdRho, dFdT, d2FdRho2, d2FdT2, d2FdRhodT, d3FdRho2dT, d3FdT2dRho, d4FdT2dRho2 ;
    
    OutFile << std::setprecision(12);

    for(double logT = logT_lo; logT<=logT_hi; logT+=dT){
        T = pow(10,logT);
        
        for(double logrho = logrho_lo; logrho<=logrho_hi; logrho+=drho){
            
            Rho = pow(10,logrho);
            
            pairs.Evaluate(Rho,T,1.);
            pairs.Evaluate_FreeEnergyDerivatives();
            
            F = pairs.u - T*pairs.s;
            dFdRho = pairs.p/Q(Rho);
            dFdT = -pairs.s;
            d2FdRho2 = pairs.dpdrho/pow(Rho,2) - 2*pairs.p/pow(Rho,3);
            d2FdT2 = -pairs.dsdT;
            d2FdRhodT = pairs.dpdT/Q(Rho);
            d3FdRho2dT = pairs.d3fdrho2dT;
            d3FdT2dRho = pairs.d3fdT2drho;
            d4FdT2dRho2 = pairs.d4fdrho2dT2;

            
            OutFile << logrho << tab << logT << tab
            << F << tab
            << dFdRho << tab
            << dFdT << tab
            << d2FdRho2 << tab
            << d2FdT2 << tab
            << d2FdRhodT << tab
            << d3FdRho2dT  << tab
            << d3FdT2dRho << tab
            << d4FdT2dRho2  << std::endl;
        }
        
        OutFile << std::endl;
    }  
    
} 
#endif


#if test1b // Phase diagram
{
    double logrho_lo = -6, logrho_hi = 11;
    double drho = 0.1;
    
    string OutName = "output/eos/phasediagram.dat";
    ofstream OutFile(OutName.c_str());
    
    OutFile << "# Phase diagram for carbon compostion "<< endl;
    OutFile << "# Rho   T(Gamma=1)  T (Gamma=175)    T(theta=1) " << std::endl;
    
    FermiDirac FD;
    
    double Rho, n, a, Omegap;

    using PhysConstants::four_pi_o3;
    using PhysConstants::four_pi;
    using PhysConstants::N_A;
    using PhysConstants::k_B;
    using PhysConstants::M_e;
    using PhysConstants::c;
    using PhysConstants::hbar;
    using PhysConstants::h;
    using PhysConstants::pi;

    double e = 4.8033e-10;

    for(double logrho = logrho_lo; logrho<=logrho_hi; logrho+=drho){
            
      Rho = pow(10,logrho);
            
      n = Rho*N_A/12;

      a = pow(four_pi_o3*n,-1/3.);

      Omegap = sqrt(four_pi*n*Q(6)*Q(e)*N_A/12);

      OutFile << logrho << tab << log10(Q(6*e)/(a*k_B*1)) << tab
            << log10(Q(6*e)/(a*k_B*175)) << tab << log10(hbar*Omegap/k_B)
            << std::endl;
       
    }  

    OutFile.close();

    OutFile.open("output/eos/degenerateline.dat");

    double T, beta, ne;

    double B_mu = 8.*pi*C(M_e)*C(c)/(3.*N_A*C(h));

    for(double logT = 3; logT<=11; logT+=0.1){
            
      T = pow(10,logT);

      beta = k_B*T/(M_e*Q(c));

      ne = 3*sqrt(2)*N_A*B_mu*pow(beta,1.5)*(FD.F(0.5,-2.5,beta) + beta*FD.F(1.5,-2.5,beta));

      OutFile << log10(2*ne/N_A) << tab << logT
            << std::endl;
       
    }  

    OutFile.close();

    OutFile.open("output/eos/degeneraterelativisticline.dat");

    ne = (1/(3*Q(pi)))*pow(M_e*c/hbar,3);

    for(double logT = 3; logT<=log10(M_e*Q(c)/k_B); logT+=0.01){

      OutFile << log10(2*ne/N_A) << tab << logT
            << std::endl;
       
    }

    OutFile.close();

    OutFile.open("output/eos/nondegeneraterelativisticline.dat");

    for(double logrho = -5; logrho<=log10(2*ne/N_A); logrho+=0.01){

      OutFile << logrho << tab << log10(M_e*Q(c)/k_B)
            << std::endl;
       
    }
    
} 
#endif

#if test2 // Calcula as qtds termodinâmicas para diferentes temperaturas
{
    std::cerr << "Test 2: Teste para qtds termodinâmicas para diferentes temperaturas" << std::endl;
    
        double T[8] = {1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11};
        
        string OutName = "output/test/eos/pressure-electron.dat";
        ofstream OutFile(OutName.c_str());
        
        string OutName2 = ".output/test/eos/energydensity-electron.dat";
        ofstream OutFile2(OutName2.c_str());
        
        string OutName3 = "output/test/eos/entropydensity-electron.dat";
        ofstream OutFile3(OutName3.c_str());
        
        string OutName4 = "output/test/eos/dPdrho-electron.dat";
        ofstream OutFile4(OutName4.c_str());
        
        string OutName5 = "output/test/eos/dUdrho-electron.dat";
        ofstream OutFile5(OutName5.c_str());
        
        string OutName6 = "output/test/eos/dPdT-electron.dat";
        ofstream OutFile6(OutName6.c_str());
        
        string OutName7 = "output/test/eos/dUdT-electron.dat";
        ofstream OutFile7(OutName7.c_str());
        
        string OutName8 = "output/test/eos/relation1-electron.dat";
        ofstream OutFile8(OutName8.c_str());
        
        string OutName9 = "output/test/eos/relation2-electron.dat";
        ofstream OutFile9(OutName9.c_str());
        
        string OutName10 = "output/test/eos/relation3-electron.dat";
        ofstream OutFile10(OutName10.c_str());
    
        string OutName11 = "output/test/eos/soundspeed-electron.dat";
        ofstream OutFile11(OutName11.c_str());
    
        OutFile << "#rho  P" << std::endl;
        OutFile4 << "#rho  dPdrho" << std::endl;
        OutFile5 << "#rho   dUdrho" << std::endl;
        OutFile6 << "#rho   dPdT" << std::endl;
        OutFile7 << "#rho   dUdT" << std::endl;
        OutFile8 << "#rho   P - rho^2 dEdrho - T dPdT" << std::endl;
        OutFile9 << "#rho   dUdT - T dSdT" << std::endl;
        OutFile10 << "#rho   dSdrho + rho^-2 PdT" << std::endl;
        OutFile11 << "#rho   cs(cm/s)" << std::endl;
        
        double rho;
    
        EoS eos;
    
        std::vector<double> X(7);
        X[1] = 1.0;
        
        for(double logrho = -4; logrho < 11; logrho+=0.1){
            
            rho = pow(10,logrho);
            
            OutFile << rho << " ";
            OutFile2 << pow(10,logrho) << " ";
            OutFile3 << pow(10,logrho) << " ";
            OutFile4 << pow(10,logrho) << " ";
            OutFile5 << pow(10,logrho) << " ";
            OutFile6 << pow(10,logrho) << " ";
            OutFile7 << pow(10,logrho) << " ";
            OutFile8 << pow(10,logrho) << " ";
            OutFile9 << pow(10,logrho) << " ";
            OutFile10 << pow(10,logrho) << " ";
            OutFile11 << pow(10,logrho) << " ";
            
            for(unsigned int i = 0; i < 8; i++){
                
                eos.Evaluate(rho,T[i],Y);
                
                OutFile << eos.P << " ";
                OutFile2 << eos.u << " ";
                OutFile3 << eos.s << " ";
                OutFile4 << eos.dPdrho << " ";
                OutFile5 << eos.dudrho << " ";
                OutFile6 << eos.dPdT << " ";
                OutFile7 << eos.dudT << " ";
                OutFile8 << fabs(eos.P - Q(rho)*eos.dudrho - T[i]*eos.dPdT)/eos.P  << " ";
                OutFile9 << fabs(eos.dudT - T[i]*eos.dsdT)/eos.dudT  << " ";
                OutFile10 << fabs(eos.dsdrho + eos.dPdT/Q(rho))/fabs(eos.dPdT/Q(rho))  << " ";
                OutFile11 << eos.cs  << " ";
            }
            
            OutFile << endl;
            OutFile2 << endl;
            OutFile3 << endl;
            OutFile4 << endl;
            OutFile5 << endl;
            OutFile6 << endl;
            OutFile7 << endl;
            OutFile8 << endl;
            OutFile9 << endl;
            OutFile10 << endl;
            OutFile11 << endl;
        }
    
    
    eos.Evaluate(2.03836e+06,1.44105e+06,Y);
    
    std::cerr << eos.dudT << std::endl;
    
        OutFile.close(); OutFile2.close(); OutFile3.close(); OutFile4.close();
        OutFile5.close(); OutFile6.close(); OutFile7.close();OutFile8.close();
        OutFile9.close(); OutFile10.close(); OutFile11.close();
} 
#endif
    
#if test3 //Teste para Nuclear Network iso7
    {
        std::cerr << "Test 3: Teste para Nuclear Network iso7" << std::endl;
        
        using PhysConstants::yr;
        
        const double tol = 1.0e-8, rtol=1.e-6, hmin=0.;
        
        std::vector<double> X(7), Xi;
        
        double h0=1e-9;
        
        Network network;
        
        X[0] = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = 0.;
        X[1] = 0.5;
        X[2] = 0.5;
        
        std::cerr << "==============================" << std::endl;
        std::cerr << "Carbon-Oxygen Explosive Burning" << std::endl;
        std::cerr << "rho = 1e9 g/cm3 and T = 3e9 K" << std::endl;
        std::cerr << "==============================" << std::endl;
        
        std::string FileName = "output/test/thermonuclear/carbonoxygen-explosive.dat";
        ofstream OutFile(FileName.c_str());
        
        double tnext, h = h0, t = 0.;
        double thd, enuc = 0.0;
        
        OutFile << t << " " << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << " " << X[4]<< " " << X[5] << " " << X[6] << " " << enuc << std::endl;

        double Rho = 1.0e9;
        double T = 3.0e9;

        double tf = 0.1;

        while (t<tf){
            
            tnext = t + h;
            
            thd = 446/sqrt(Rho);
            
            network.Rho = Rho*exp(-t/thd);
            network.T = T*exp(-t/(3*thd));
            
            Xi = X;
            
            Odeint<StepperSie<Network> > ode(X,t,tnext,tol,rtol,h,hmin,network);
            
            ode.integrate();

            enuc = network.EnergyGeneratedIntegrated(X,Xi,h);
            
            h = ode.s.hnext; //Atualiza o tamanho do passo

            t = tnext;
            
            OutFile << t << " " << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << " " << X[4]<< " " << X[5] << " " << X[6] << " " << enuc << std::endl;
        }
        
        OutFile.close();
        
        std::cerr << "==============================" << std::endl;
        std::cerr << "Helium Explosive Burning" << std::endl;
        std::cerr << "rho = 1e8 g/cm3 and T = 3e9 K" << std::endl;
        std::cerr << "==============================" << std::endl;
        
        FileName = "output/test/thermonuclear/helium-explosive.dat";
        OutFile.open(FileName.c_str());
        
        X[0] = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = 0.;
        X[0] = 1.0;
        
        Rho = 1.0e8;
        T = 3.0e9;
        
        t = 0.0;
        h = h0;
        
        while (t<tf){
            
            tnext = t + h;
            
            thd = 446/sqrt(Rho);
            
            network.Rho = Rho*exp(-t/thd);
            network.T = T*exp(-t/(3*thd));
            
            Xi = X;
            
            Odeint<StepperSie<Network> > ode(X,t,tnext,tol,rtol,h,hmin,network);
            
            ode.integrate();
            
            enuc = network.EnergyGeneratedIntegrated(X,Xi,h);
            
            h = ode.s.hnext; //Atualiza o tamanho do passo
            
            t = tnext;
            
            OutFile << t << " " << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << " " << X[4]<< " " << X[5] << " " << X[6] << " " << enuc << std::endl;
        }
        
        OutFile.close();
        
        std::cerr << "==============================" << std::endl;
        std::cerr << "Silicon Explosive Burning" << std::endl;
        std::cerr << "rho = 1e9 g/cm3 and T = 5e9 K" << std::endl;
        std::cerr << "==============================" << std::endl;
        
        FileName = "output/test/thermonuclear/silicon-explosive.dat";
        OutFile.open(FileName.c_str());
        
        X[0] = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = 0.;
        X[5] = 1.0;
        
        Rho = 1.0e9;
        T = 5.0e9;
        
        t = 0.0;
        h = h0;
        
        while (t<tf){
            
            tnext = t + h;
            
            thd = 446/sqrt(Rho);
            
            network.Rho = Rho*exp(-t/thd);
            network.T = T*exp(-t/(3*thd));
            
            Xi = X;
            
            Odeint<StepperSie<Network> > ode(X,t,tnext,tol,rtol,h,hmin,network);
            
            ode.integrate();
            
            enuc = network.EnergyGeneratedIntegrated(X,Xi,h);
            
            h = ode.s.hnext; //Atualiza o tamanho do passo
            
            t = tnext;
            
            OutFile << t << " " << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << " " << X[4]<< " " << X[5] << " " << X[6] << " " << enuc << std::endl;
        }
        
        
        Xi.clear();

        OutFile.close();
        
        std::cerr << "==============================" << std::endl;
        std::cerr << "NuclearStatisticalEquilibrium" << std::endl;
        std::cerr << " T = 6e9 K" << std::endl;
        std::cerr << "==============================" << std::endl;
        
        FileName = "output/test/thermonuclear/nuclearstatisticalequilibrium-T6e9.dat";
        OutFile.open(FileName.c_str());
        OutFile << "# logrho   XHe    XSi    XNi " << enuc << std::endl;

        for (double logrho = 6.0; logrho<10.0; logrho+=0.1){

            network.Nuclear_Network::NuclearStatisticalEquilibrium(pow(10.0,logrho),6e9,X);

            OutFile << logrho << " " << X[0] << " " << X[5] << " " << X[6] << " " << enuc << std::endl;
        }

        X.clear();
        
    }
#endif

#if test4 //Teste para taxa de reação de He4
{
    std::cerr << "Test 4: Teste para taxa de reação de He4" << std::endl;
    
    using PhysConstants::yr;

    const double rtol=1.e-8;
    
    std::vector<double> Y(7);
    
    Y[0] = Y[2] = Y[3] = Y[4] = Y[5] = Y[6] = Y[5] = 0.;
    Y[1] = 1.0/12;
    
    std::string FileName = "output/thermonuclear/reactionrateC.dat";
    ofstream OutFile(FileName.c_str());

    double logT1, logT2;

    rootHe func;
    func.Y = Y;

    for (double logrho = 5.5; logrho < 10; logrho+=0.1){
    
        func.Rho = pow(10,logrho);

        func.tau = 1e6*yr;

        logT1 = zbrent(func,8.0,11.0,rtol);

        func.tau = 446/sqrt(func.Rho);

        logT2 = zbrent(func,8.0,11.0,rtol);
        
        OutFile << logrho << " " << logT1 << " " << logT2 << std::endl;
        
    }
    
    Y.clear();
    
}
#endif

#if test5 //Teste para condutividade térmica por elétrons
{
    std::cerr << "Test 5: Teste para condutividade térmica por elétrons" << std::endl;
    
    ThermalConductivity cond;

    std::vector<double> Y(7);
    
    Y[0] = Y[2] = Y[3] = Y[4] = Y[5] = Y[6] = Y[5] = 0.;
    Y[1] = 1.0;

    std::string FileName = "output/test/opacity/thermalconductivity-carbon.dat";
    ofstream OutFile(FileName.c_str());

    for (double logrho = -5.5; logrho < 12.0; logrho+=0.05){
        OutFile << logrho << " " << cond.Evaluate_ThermalConductivity(logrho,5.0,Y) << " " << cond.Evaluate_ThermalConductivity(logrho,6.0,Y) << " " << cond.Evaluate_ThermalConductivity(logrho,7.0,Y) << " " << cond.Evaluate_ThermalConductivity(logrho,8.0,Y) << " " << cond.Evaluate_ThermalConductivity(logrho,9.0,Y) << " " << cond.Evaluate_ThermalConductivity(logrho,10.0,Y) << std::endl;
    }

    // std::string FileName2 = "output/test/opacity/thermalconductivity-carbon-degeneratecase.dat";
    // ofstream OutFile2(FileName2.c_str());

     using PhysConstants::hbar;
     using PhysConstants::pi;
     using PhysConstants::N_A;
     using PhysConstants::c;
     using PhysConstants::e_c;
     using PhysConstants::k_B;
     using PhysConstants::M_e;
     using PhysConstants::a;

    // for (double logT = 4; logT < 7.5; logT+=0.01){
    //     OutFile2 << logT << " " << cond.Degenerateconductivity(4.0,logT,Y) << " " << cond.Evaluate_ThermalConductivity(4.0,logT,Y) << std::endl;
    // }
    
     double kapparad;
    
     double sigmacond, kappacond, kappa, T;

     std::string FileName3 = "output/test/opacity/opacity-carbon-logT8.dat";
     ofstream OutFile3(FileName3.c_str());

     for (double logrho = -5.5; logrho < 12.0; logrho+=0.05){

         T = pow(10,8.0);

         sigmacond = pow(10,cond.Evaluate_ThermalConductivity(logrho,log10(T),Y));
         kappacond = 4*a*c*C(T)/(3*pow(10,logrho)*sigmacond);

         kapparad =  1.e23*pow(10,logrho)*pow(T,-3.5) + 0.2;

         kappa = 1./(1./kapparad+1./kappacond);

         OutFile3 << logrho << " " << kapparad << " " << kappacond << " " << kappa << std::endl;
     }
    
}
#endif


#if test10 //Test to evolve white dwarf during 1e8 years
{
  std::cerr << "Test 10: Test to evolve white dwarf during 1e8 years" << std::endl;

  omp_set_num_threads(2); 
    
  using PhysConstants::M_Sun;
  using PhysConstants::R_Sun;
  using PhysConstants::yr;

  std::vector<double> X(7);
  X[1] = 0.5;
  X[2] = 0.5;

  Star star(108,1.08*M_Sun,5e7,pow(10.0,8.),X);

  std::cerr << "The first hydrostatic equilibrium...";

  star.HydrostaticEquilibrium();

  std::cerr << "done!"<< std::endl;

  //Evolução térmica da estrela
  const double tol = 1.e-5, rtol=1.e-4, hmin=0.;

  star.OutputFilePath = "output/equilibrium/";

  //OutDynamics myout(star);

  //myout.AllVariables(0.);

  double h = yr, tnext, t = 0.;

  while (t<1e8*yr){
    tnext = t + h;

    star.Time_Integration(t,h,rtol,tol,2);
      
    star.HydrostaticEquilibrium();
  	
    t = tnext;

    std::cerr << "t = " << tnext/yr << " yr" << std::endl;

    //myout.AllVariables(tnext);
  }

  //Equilibrio hidrostatico final
  star.HydrostaticEquilibrium();

  double E = star.Gravitational_Energy() + star.Internal_Energy();
    
    std::string FileName = star.OutputFilePath + "white-dwarf-M" + to_string(star.Mass/M_Sun) + "-logTc" + to_string(log10(star.shell[0].T) ) + "-N" + to_string(star.Nshells) + ".dat";
    ofstream OutFile(FileName.c_str());
    OutFile << "# M    R    Rdot     T      Y_He     Y_C    Y_O    Y_Ne    Y_Mg    Y_Si    Y_Ni" << std::endl;
    OutFile << "# Energy = " << E << " ergs" << std::endl;
    
    for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
        OutFile << star.shell[j].M << " " << star.shell[j].Rout << " " << 0. << " " << star.shell[j].T << " " << star.shell[j].Y[0] << " " << star.shell[j].Y[1] << " " << star.shell[j].Y[2] << " " << star.shell[j].Y[3] << " " << star.shell[j].Y[4] << " " << star.shell[j].Y[5] << " " << star.shell[j].Y[6] << std::endl;
    }

    std::string FileName2 = star.OutputFilePath + "initial-condition-M" + to_string(star.Mass/M_Sun) + "-logTc" + to_string(log10(star.shell[0].T) ) + "-N" + to_string(star.Nshells) + ".dat";
    ofstream OutFile2(FileName2.c_str());
    OutFile2 << "# R    Rho       T" << std::endl;
    OutFile2 << "# Energy = " << E << " ergs" << std::endl;
    
    for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
        OutFile2 << star.shell[j].Rout << " " << star.shell[j].Rho << " " << star.shell[j].T << std::endl;
    }

}
#endif


#if test12 //Free-fall test
{
  std::cerr << "Test 12: Free-fall test" << std::endl;
    
    using PhysConstants::G;
    using PhysConstants::pi;
    using PhysConstants::four_pi_o3;
    using PhysConstants::yr;

    std::vector<double> X(7);
    X[1] = 0.5;
    X[2] = 0.5;

    double R0 = 1e9;
    double rho0 = 1e7;

    double M = rho0*four_pi_o3*C(R0);
    
    Star star(100,M,rho0,1e7,X);

    star.eos.id_eos = 0;
    star.art_viscosity = true;
    star.reactions = false;
    star.luminosity = false;

    star.Update_Interior();
    
    star.OutputFilePath = "output/test/";
    //star.InputFilePath = "../input/";
  
    const double tol = 1.e-6, rtol=1.e-4;
  
    double tnext, dt = 0.001, t = 0., tlast=t;

    std::string FileName = star.OutputFilePath + "freefall.dat";
    ofstream OutFile(FileName.c_str());

    std::string FileName2 = star.OutputFilePath + "freefall-rho.dat";
    ofstream OutFile2(FileName2.c_str());


    std::string FileName3 = star.OutputFilePath + "freefall-profile-t=0.dat";
    ofstream OutFile3(FileName3.c_str());

    for (unsigned int j = 0; j < star.Nshells; j++){
        OutFile3 << star.shell[j].Rin << " " << star.shell[j].Rho << std::endl;
        OutFile3 << star.shell[j].Rout << " " << star.shell[j].Rho << std::endl;
    }
    OutFile3 << star.shell[star.Nshells-1].Rout << " " << 0. << std::endl;

    OutFile << star.shell[star.Nshells-1].Rout/R0 << " " << pi/2.-t*sqrt(8*pi*G*rho0/3.) << " " << pi/2. -(sqrt(1-star.shell[star.Nshells-1].Rout/R0)*sqrt(star.shell[star.Nshells-1].Rout/R0) + asin(sqrt(1-star.shell[star.Nshells-1].Rout/R0)) ) << std::endl;

    OutFile2 << t << " " << star.shell[star.Nshells-1].Rho << " " << star.shell[star.Nshells-1].Rout << std::endl;
    
    double ti[5] = {0.550,0.645,0.660,0.6635,0.6641};
    unsigned int i = 0;

    unsigned int tictac = 0;

    while (t < 0.66431){
      
        tnext = t + dt;

        star.Time_Integration(t,dt,rtol,tol,1);

        tictac++;

        if (tictac == 50){
            tictac = 0;
        
            std::cerr << "t =" << t << " s" << std::endl;

            OutFile << star.shell[star.Nshells-1].Rout/R0 << " " << pi/2.-t*sqrt(8*pi*G*rho0/3.) << " " << pi/2. -(sqrt(1-star.shell[star.Nshells-1].Rout/R0)*sqrt(star.shell[star.Nshells-1].Rout/R0) + asin(sqrt(1-star.shell[star.Nshells-1].Rout/R0)) ) << std::endl;

            OutFile2 << t << " " << star.shell[star.Nshells-1].Rho << " " << star.shell[star.Nshells-1].Rout << std::endl;
        }
        
        
        if (t > ti[i]){
            string FileNamex = star.OutputFilePath + "freefall-profile-t=" + to_string(t) + ".dat";
            ofstream OutFilex(FileNamex.c_str());
            
            for (unsigned int j = 0; j < star.Nshells; j++){
                OutFilex << star.shell[j].Rin << " " << star.shell[j].Rho << std::endl;
                OutFilex << star.shell[j].Rout << " " << star.shell[j].Rho << std::endl;
            }
            OutFilex << star.shell[star.Nshells-1].Rout << " " << 0. << std::endl;
            i++;
        }
    }

    FileName3 = star.OutputFilePath + "freefall-profile-t=" + to_string(t) + ".dat";
    ofstream OutFile4(FileName3.c_str());

    for (unsigned int j = 0; j < star.Nshells; j++){
        OutFile4 << star.shell[j].Rin << " " << star.shell[j].Rho << std::endl;
        OutFile4 << star.shell[j].Rout << " " << star.shell[j].Rho << std::endl;
    }
    OutFile4 << star.shell[star.Nshells-1].Rout << " " << 0. << std::endl;
    
    std::cerr << "done!" << std::endl;

}
#endif

#if test13 //Sod's Shock Wave test
{
  std::cerr << "Test 13: Sod's Shock Wave test" << std::endl;
    
    using PhysConstants::G;
    using PhysConstants::pi;
    using PhysConstants::four_pi_o3;
    using PhysConstants::M_Sun;

    std::vector<double> X(7);
    X[1] = 0.5;
    X[2] = 0.5;

    double rho0 = 1e7;
    
    Star star(1000,M_Sun,rho0,1e4,X);

    star.eos.id_eos = 2;
    star.art_viscosity = true;
    star.reactions = false;
    star.luminosity = false;
    star.gravity = false;
    star.tau = 0.0;

    star.shell[0].Rin = 0;
    star.shell[0].Rout = pow(star.shell[0].M/(four_pi_o3*rho0),1/3.);
    for (unsigned int j = 1; j < star.Nshells/5; j++){
        star.shell[j].Rin = star.shell[j-1].Rout;
        star.shell[j].Rout = pow(C(star.shell[j].Rin) + star.shell[j].M/(four_pi_o3*rho0),1/3.);
    }

    for (unsigned int j = star.Nshells/5; j < star.Nshells; j++){
        star.shell[j].Rin = star.shell[j-1].Rout;
        star.shell[j].Rout = pow(C(star.shell[j].Rin) + star.shell[j].M/(four_pi_o3*0.125*rho0),1/3.);
    }

    star.Update_Interior();
    
    star.OutputFilePath = "output/test/sodshock/";
  
    const double tol = 1.e-6, rtol=1.e-4;
  
    double tnext, dt = 0.001, t = 0., tlast=t;

    std::string FileName = star.OutputFilePath + "sod-external.dat";
    ofstream OutFile(FileName.c_str());

    std::string FileName1 = star.OutputFilePath + "sod-internal.dat";
    ofstream OutFile1(FileName1.c_str());

    OutFile << t << " " << star.shell[star.Nshells-1].Rout << " " << star.shell[star.Nshells-1].Rho << " " << star.shell[star.Nshells-1].P << std::endl;

    OutFile1 << t << " " << star.shell[0].Rout << " " << star.shell[0].Rho << " " << star.shell[0].P << std::endl;

    std::string FileName3 = star.OutputFilePath + "sod-profile-t=0.dat";
    ofstream OutFile3(FileName3.c_str());

    OutFile3 << "# R(cm)   Rho (g/cm3)   P (erg/cm3)   vR (cm/s) " << std::endl;

    OutFile3 << 0 << " " << star.shell[0].Rho << " " << star.shell[0].P << " " << star.shell[0].Rdot_out << std::endl;
    for (unsigned int j = 0; j < star.Nshells; j++){
        OutFile3 << star.shell[j].Rout << " " << star.shell[j].Rho  << " " << star.shell[j].P << " " << star.shell[j].Rdot_out << std::endl;
    }
    OutFile3 << star.shell[star.Nshells-1].Rout << " " << 0. << " " << 0 << "  " << 0 <<  std::endl;

    OutFile3.close();

    bool print = true;

    while (t < 0.3){
      
        tnext = t + dt;

        star.Time_Integration(t,dt,rtol,tol,1);

        std::cerr << "t =" << t << " s" << std::endl;

        if (t - tlast > 0.001){ 
            tlast = t;       
            std::cerr << "t =" << t << " s" << std::endl;

            OutFile << t << " " << star.shell[star.Nshells-1].Rout << " " << star.shell[star.Nshells-1].Rho << " " << star.shell[star.Nshells-1].P << std::endl;

            OutFile1 << t << " " << star.shell[0].Rout << " " << star.shell[0].Rho << " " << star.shell[0].P << std::endl;
        }

        if (t > .15 && print == true){
            print = false;
            FileName3 = star.OutputFilePath + "sod-profile-medium.dat";
            OutFile3.open(FileName3.c_str());

            OutFile3 << 0 << " " << star.shell[0].Rho << " " << star.shell[0].P << " " << star.shell[0].Rdot_out << std::endl;
            for (unsigned int j = 0; j < star.Nshells; j++){
                OutFile3 << star.shell[j].Rout << " " << star.shell[j].Rho  << " " << star.shell[j].P << " " << star.shell[j].Rdot_out << std::endl;
            }
            OutFile3 << star.shell[star.Nshells-1].Rout << " " << 0. << " " << 0 << "  " << 0 <<  std::endl;
            OutFile3.close();
        }
    }
    
    std::cerr << "done!" << std::endl;

    FileName3 = star.OutputFilePath + "sod-profile-final.dat";
    OutFile3.open(FileName3.c_str());

    OutFile3 << 0 << " " << star.shell[0].Rho << " " << star.shell[0].P << " " << star.shell[0].Rdot_out << std::endl;
    for (unsigned int j = 0; j < star.Nshells; j++){
        OutFile3 << star.shell[j].Rout << " " << star.shell[j].Rho  << " " << star.shell[j].P << " " << star.shell[j].Rdot_out << std::endl;
    }
    OutFile3 << star.shell[star.Nshells-1].Rout << " " << 0. << " " << 0 << "  " << 0 <<  std::endl;

}
#endif

#if test14 //Sedov's Blast Wave test
{
  std::cerr << "Test 14: Sedov's Blast Wave test" << std::endl;
    
    using PhysConstants::G;
    using PhysConstants::pi;
    using PhysConstants::four_pi_o3;
    using PhysConstants::yr;
    using PhysConstants::M_Sun;

    std::vector<double> X(7);
    X[1] = 0.5;
    X[2] = 0.5;

    double rho0 = 1.0e7;
    
    Star star(1000,M_Sun,rho0,1e4,X);

    star.eos.id_eos = 2;
    star.art_viscosity = true;
    star.reactions = false;
    star.luminosity = false;
    star.gravity = false;

    //Increase the temperature in the center
    star.Update_Interior();

    star.shell[0].Eps = 10000*star.shell[0].u;
    
    star.OutputFilePath = "output/test/sedov/";
  
    const double tol = 1.e-6, rtol=1.e-4;
  
    double tnext, dt = 1.e-5, t = 0., tlast=t;

    std::string FileName3 = star.OutputFilePath + "sedov-profile-initial.dat";
    ofstream OutFile3(FileName3.c_str());

    OutFile3 << 0 << " " << star.shell[0].Rho << " " << star.shell[0].P << " " << star.shell[0].Rdot_out << std::endl;
    for (unsigned int j = 0; j < star.Nshells; j++){
        OutFile3 << star.shell[j].Rout << " " << star.shell[j].Rho  << " " << star.shell[j].P << " " << star.shell[j].Rdot_out << std::endl;
    }
    OutFile3 << star.shell[star.Nshells-1].Rout << " " << 0. << " " << 0 << "  " << 0 <<  std::endl;

    OutFile3.close();

    bool print = true;

    while (t < 0.002){
      
        tnext = t + dt;

        star.Time_Integration(t,dt,rtol,tol,3);
          
        t = tnext;

        if (t - tlast > 0.00){ 
            tlast = t;       
            std::cerr << "t =" << t << " s" << std::endl;
        }

        if (t > 0.001 && print == true){
            print = false;
            FileName3 = star.OutputFilePath + "sedov-profile-medium.dat";
            OutFile3.open(FileName3.c_str());

            OutFile3 << 0 << " " << star.shell[0].Rho << " " << star.shell[0].P << " " << star.shell[0].Rdot_out << std::endl;
            for (unsigned int j = 0; j < star.Nshells; j++){
                OutFile3 << star.shell[j].Rout << " " << star.shell[j].Rho  << " " << star.shell[j].P << " " << star.shell[j].Rdot_out << std::endl;
            }
            OutFile3 << star.shell[star.Nshells-1].Rout << " " << 0. << " " << 0 << "  " << 0 <<  std::endl;
            OutFile3.close();
        }
    }
    
    std::cerr << "done!" << std::endl;

    FileName3 = star.OutputFilePath + "sedov-profile-final.dat";
    OutFile3.open(FileName3.c_str());

    OutFile3 << 0 << " " << star.shell[0].Rho << " " << star.shell[0].P << " " << star.shell[0].Rdot_out << std::endl;
    for (unsigned int j = 0; j < star.Nshells; j++){
        OutFile3 << star.shell[j].Rout << " " << star.shell[j].Rho  << " " << star.shell[j].P << " " << star.shell[j].Rdot_out << std::endl;
    }
    OutFile3 << star.shell[star.Nshells-1].Rout << " " << 0. << " " << 0 << "  " << 0 <<  std::endl;

}
#endif


#if test20 //Dynamics with time delay given by tau
{
  std::cerr << "Test 20: Dynamics with time delay given by tau" << std::endl;
    
    using PhysConstants::M_Sun;
    using PhysConstants::R_Sun;
    using PhysConstants::G;
    using PhysConstants::four_pi;
    using PhysConstants::four_pi_o3;
    using PhysConstants::yr;

    Star star;
    star.OutputFilePath = "output/dynamics/";
    star.tau = 0;
  
    const double tol = 1.e-6, rtol=1.e-4;
  
    double tnext, dt, t = 0., tlast=0.; bool stop = false;
    
    double Mdot = 1.e-6*M_Sun/yr;
    
      star.shell.resize(star.Nshells+1);
      star.Nshells = star.shell.size();
    
      star.shell[star.Nshells-1].M = 0.001*M_Sun;
    
     star.shell[star.Nshells-1].Mdot = Mdot;
    
      star.shell[star.Nshells-1].Y[0] = 1.0/4.;
      star.shell[star.Nshells-1].Yi = star.shell[star.Nshells-1].Y;
      star.shell[star.Nshells-1].T = star.shell[star.Nshells-2].T;
    
      star.shell[star.Nshells-1].Rin = star.shell[star.Nshells-2].Rout;
      star.shell[star.Nshells-1].Rout = pow(C(star.shell[star.Nshells-1].Rin) + star.shell[star.Nshells-1].M/(four_pi_o3*star.shell[star.Nshells-2].Rho),1/3.);
    
     star.shell[star.Nshells-1].Rdot_in = 0.;
     star.shell[star.Nshells-1].Rdot_out = -Mdot/(four_pi*Q(star.shell[star.Nshells-2].Rho));

    star.Update_Interior();
    
    dt = 10*yr;

    std::vector<double> Rmem(200);

    std::string FileName0 = star.OutputFilePath + "shell-accreted.dat";
    ofstream OutFile0(FileName0.c_str());

    std::string FileName2 = star.OutputFilePath + "shell-central.dat";
    ofstream OutFile2(FileName2.c_str());
    
    std::string FileName = star.OutputFilePath + "whitedwarf-accreted.dat";

    OutFile0 << 0 << " " << log10(star.shell[107].Rho) << " " << log10(star.shell[107].T) << std::endl;
    OutFile2 << 0 << " " << log10(star.shell[0].Rho) << " " << log10(star.shell[0].T) << std::endl;
    
    rootHe func;

    while (t < 1e8*yr){
      
        tnext = t + dt;

        star.shell[star.Nshells-1].M += Mdot*dt;

        if (star.shell[star.Nshells-1].M > 0.012*M_Sun) {

            std::cerr << "======================= " << std::endl;
            std::cerr << "Number of shells = " << star.Nshells+1 << std::endl;

            star.shell.resize(star.Nshells+1);
            star.Nshells = star.shell.size();

            star.shell[star.Nshells-1].M = star.shell[star.Nshells-2].M - 0.01*M_Sun;
            star.shell[star.Nshells-2].M = 0.01*M_Sun;

            star.shell[star.Nshells-1].Mdot = Mdot;
            star.shell[star.Nshells-2].Mdot = 0.0;

            star.shell[star.Nshells-1].Y[0] = 1.0/4.;
            star.shell[star.Nshells-1].Yi = star.shell[star.Nshells-1].Y;
            star.shell[star.Nshells-1].T = star.shell[star.Nshells-2].T;
        
            star.shell[star.Nshells-1].Rin = star.shell[star.Nshells-2].Rout;
            star.shell[star.Nshells-1].Rout = pow(C(star.shell[star.Nshells-1].Rin) + star.shell[star.Nshells-1].M/(four_pi_o3*star.shell[star.Nshells-2].Rho),1/3.);
        
            star.shell[star.Nshells-1].Rdot_in = star.shell[star.Nshells-2].Rdot_out;
            star.shell[star.Nshells-1].Rdot_out = -Mdot/(four_pi*Q(star.shell[star.Nshells-2].Rho));
            
            star.Update_Interior();
            
            std::string FileName = star.OutputFilePath + "wdprofile-t" + to_string(t/yr) + "-Mdot" + to_string(Mdot/(M_Sun/yr) ) + ".dat";
            ofstream OutFile(FileName.c_str());
            OutFile << "# R    Rho     T" << std::endl;
            OutFile << "# Energy = " << star.Gravitational_Energy() + star.Internal_Energy() << " ergs" << std::endl;
            
            for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
                OutFile << star.shell[j].Rout << " " << star.shell[j].Rho << " " << star.shell[j].T << std::endl;
            }
            
            std::string FileName1 = star.OutputFilePath + "whitedwarf-t" + to_string(t/yr) + "-Mdot" + to_string(Mdot/(M_Sun/yr) ) + ".dat";
            ofstream OutFile1(FileName1.c_str());
            OutFile1 << "# M    R    Rdot     T      Y_He     Y_C    Y_O    Y_Ne    Y_Mg    Y_Si    Y_Ni" << std::endl;
            OutFile1 << "# Energy = " << star.Gravitational_Energy() + star.Internal_Energy() << " ergs" << std::endl;
            
            for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
                OutFile1 << star.shell[j].M << " " << star.shell[j].Rout << " " << star.shell[j].Rdot_out << " " << star.shell[j].T << " " << star.shell[j].Y[0] << " " << star.shell[j].Y[1] << " " << star.shell[j].Y[2] << " " << star.shell[j].Y[3] << " " << star.shell[j].Y[4] << " " << star.shell[j].Y[5] << " " << star.shell[j].Y[6] << std::endl;
            }
        }
        
        else{ 
            #pragma omp parallel for
            for (std::vector<int>::size_type j = 0; j < star.shell.size(); j++){
                Rmem[j] = star.shell[j].Rout;
            }
            
            star.HydrostaticEquilibrium();

            #pragma omp parallel for
            for (std::vector<int>::size_type j = 0; j < star.shell.size(); j++){
                star.shell[j].Rdot_out = (star.shell[j].Rout - Rmem[j])/dt;
            }

            star.shell[0].Rdot_in = 0.;
            #pragma omp parallel for
            for (std::vector<int>::size_type j = 1; j < star.shell.size(); j++){
                star.shell[j].Rdot_in = star.shell[j-1].Rdot_out;
            }

            star.Update_Interior();
        }

        if (star.reactions) star.NuclearReactions(dt,rtol,tol);
          
        star.Time_Integration(t,dt,rtol,tol,2);

        dt = MIN(dt,0.001*M_Sun/Mdot);

        if (tnext-tlast > 0.0001*M_Sun/Mdot) {
            tlast = tnext;
            
            std::cerr << "t = " << tnext/yr << " yr" << std::endl;
            
            OutFile0 << tnext/yr << " " << log10(star.shell[107].Rho) << " " << log10(star.shell[107].T) << std::endl;
            OutFile2 << tnext/yr << " " << log10(star.shell[0].Rho) << " " << log10(star.shell[0].T) << std::endl;
            
            ofstream OutFile(FileName.c_str());
            OutFile << "# R    Rho     T" << std::endl;
            OutFile << "# Energy = " << star.Gravitational_Energy() + star.Internal_Energy() << " ergs" << std::endl;
            
            for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
                OutFile << star.shell[j].Rout << " " << star.shell[j].Rho << " " << star.shell[j].T << std::endl;
            }
            
        }
        
        if (star.Nshells > 108){
            
            for (std::vector<int>::size_type j = 108; j < star.shell.size(); j++){
                
                if (star.shell[j].T > pow(10,8.0) && star.shell[j].Rho > pow(10,5.5) ){
                    
                    func.Y = star.shell[j].Y;
                    func.Rho = star.shell[j].Rho;
                    func.tau = 3*446/sqrt(func.Rho);
                    double logT1 = zbrent(func,8.0,11.0,rtol);
                    
                    if (star.shell[j].T > pow(10,logT1) ){
                        stop = true;
                        
                        std::cerr << "tau = " << 446/sqrt(star.shell[j].Rho) << std::endl;
                    }
                }
            }

            for (std::vector<int>::size_type j = 0; j < 108; j++){
                
                if (star.shell[j].T > pow(10,9.0) && star.shell[j].Rho > pow(10,6.0) ){
                    
                    func.Y = star.shell[j].Y;
                    func.Rho = star.shell[j].Rho;
                    func.tau = 3*446/sqrt(func.Rho);
                    double logT1 = zbrent(func,9.0,10.0,rtol);
                    
                    if (star.shell[j].T > pow(10,logT1) ){
                        stop = true;
                        
                        std::cerr << "tau = " << 446/sqrt(star.shell[j].Rho) << std::endl;
                    }
                }
            }
        }

        if (stop) break; 
          
        t = tnext;

    } 

    std::string FileName3 = star.OutputFilePath + "wd-accreted-M" + to_string(star.Mass/M_Sun) + "-logTc" + to_string(log10(star.shell[0].T) ) + "-N" + to_string(star.Nshells) + ".dat";
    ofstream OutFile3(FileName3.c_str());
    OutFile3 << "# M    R    Rdot     T      Y_He     Y_C    Y_O    Y_Ne    Y_Mg    Y_Si    Y_Ni" << std::endl;
    OutFile3 << "# Energy = " << star.Gravitational_Energy() + star.Internal_Energy() << " ergs" << std::endl;
    
    for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
        OutFile3 << star.shell[j].M << " " << star.shell[j].Rout << " " << star.shell[j].Rdot_out << " " << star.shell[j].T << " " << star.shell[j].Y[0] << " " << star.shell[j].Y[1] << " " << star.shell[j].Y[2] << " " << star.shell[j].Y[3] << " " << star.shell[j].Y[4] << " " << star.shell[j].Y[5] << " " << star.shell[j].Y[6] << std::endl;
    }
    
    std::cerr << "done!" << std::endl;

}
#endif

#if test21 //Dynamics with artificial mass pertubation
{
  std::cerr << "Test 21: Dynamics with artificial mass pertubation" << std::endl;
    
    using PhysConstants::M_Sun;
    using PhysConstants::R_Sun;
    using PhysConstants::G;
    using PhysConstants::four_pi;
    using PhysConstants::four_pi_o3;
    using PhysConstants::yr;

    std::vector<double> X(7);
    X[1] = 0.5;
    X[2] = 0.5;
    
    Star star(108,1.08*M_Sun,5e7,pow(10.0,7.53),X);
    
    star.Update_Interior();
    
    star.OutputFilePath = "output/dynamics/";
    star.tau = 0.01;
    
    star.shell.resize(star.Nshells+8);
    star.Nshells = star.shell.size();
    
    for (unsigned int j = 108; j < star.Nshells; j++){
    
        star.shell[j].M = 0.001*M_Sun;
    
        star.shell[j].Mdot = 0.;
    
        star.shell[j].Y[0] = 1.0/4.;
        star.shell[j].Yi = star.shell[j].Y;
        star.shell[j].T = pow(10.,7.9);
    
        star.shell[j].Rin = star.shell[j-1].Rout;
        star.shell[j].Rout = pow(C(star.shell[j].Rin) + star.shell[j].M/(four_pi_o3*star.shell[107].Rho),1/3.);
    
        star.shell[j].Rdot_in = 0.;
        star.shell[j].Rdot_out = 0.;
    }
    
    for (unsigned int j = 0; j < 113; j++){
        star.shell[j].T = pow(10.,(7.53+(8.5-7.53)*(j/112.)));
    }
    
    for (unsigned int j = 113; j < star.Nshells; j++){
        star.shell[j].T = pow(10.,(8.5-(8.5-7.7)*(j-113)/3.));
    }
    
    star.Update_Interior();
    
    //star.HydrostaticEquilibrium();
    
    std::cerr << "Condição Inicial criada!" << std::endl;
  
    const double tol = 1.e-10, rtol=1.e-8;
  
    double tnext, dt = 1.e-6, t = 0., tlast=t;

    OutDynamics myout(star);
    myout.AllVariables(0.);

    star.reactions = true;

    while (t < 100){
      
        tnext = t + dt;
          
        star.Time_Integration(t,dt,rtol,tol,3);

        t = tnext;

        if (t - tlast > 0.00){

            tlast = t;

            myout.AllVariables(t);
          
            std::cerr << "t =" << t << " s" << std::endl;
        }
    }
    
    std::cerr << "done!" << std::endl;

}
#endif

	end = time(0);


	std::cerr << "Time elapsed: " << difftime(end,begin) << " s"<< std::endl;
    
	return 0;
}
