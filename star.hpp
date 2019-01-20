//=====================================================================
// Author: Elvis Soares
// Date: September 2013
// last modification: Sep 2013
//=====================================================================
#ifndef _STAR_HPP_
#define _STAR_HPP_

#include <vector>
#include <string>
#include <fstream>

#include "eos.hpp"
#include "cpp_constants_CGS.hpp"
#include "thermonuclear.hpp"
#include "conductivity.hpp"

double TemperaturebyEntropy(const double &s, const double &rho, const std::vector<double> &X);

struct Shells {
    
    Nuclear_Network network;
    ThermalConductivity cond;

    double M, Rout, Rin, Rdot_out, Rdot_in, T; //Quantidades definidas por Star
    std::vector<double> X, Xi;

    double Vol, Rho, Gamma, q, Q, Eps, Qd, kappa, L, Mdot; //Quantidades calculadas aqui
    double ksi, ksidot;
    double Vg, U;

    double P, u, s, Cv, Cp, cs, dUdrho; // Quantidades termodinâmcias

    double Voldot, Rhodot;

    double dMdt, dRdt, dRdotdt, dTdt, dQddt, dUdt;
    std::vector<double> dXdt;
    std::vector< std::vector<double> > dFdX;

    Shells(){
        X.resize(7);
        Xi.resize(7);
        dXdt.resize(7);
        dFdX.resize(7);
        for (std::vector<int>::size_type i = 0; i < 7; i++)
            dFdX[i].resize(7);

        M = 0.0; Rout = 0.0; Rin = 0.0; Rdot_out=0.0; Rdot_in=0.0; T=0.0; 
        dMdt=0.0; dRdt=0.0; dRdotdt=0.0; dTdt=0.0; dQddt=0.0, dUdt = 0.0;

        Voldot = 0.0; Rhodot = 0.0;
        Mdot = 0.;
        Qd = 0.;
        Eps = 0.;
        L = 0.;
        q = 0.;
        Q = 0.0;
    };

    void Set_Mass(const double &Mass);
    void Set_Radius(const double &Radius, const double &dRadius);
    void Set_RadialVelocity(const double &RadialVelocity, const double &dRadialVelocity);
    void Set_Temperature(const double &Temperature);
    void Set_Abundance(const std::vector<double> &Xin);

    void Calculate_Ksi();

    void Calculate_Volume();
    void Calculate_Density();
    
    void Calculate_ArtificialViscosity();
    void Calculate_Opacity();
};

//=====================================================================
// Here the classes for Star are defined. It consists of a
// collection of all informations about its interior
//=====================================================================

class Star {
    
public :

    EoS eos;

    static std::string ConfigFileName, OutputFilePath, InputFilePath, OutputFileName, EndFileName, InitialConditionFileName;

    bool reactions, luminosity, art_viscosity, gravity;

    unsigned int Nshells;

    double Mass, Density, Temperature, Luminosity;
    std::vector<double> X;

    double tau, dt_integrator;

    std::vector<Shells> shell;

    Star(int nshells, double mass, double density, double temp, std::vector<double> Xin): Nshells(nshells), Mass(mass), Density(density), Temperature(temp), X(Xin){ // class constructor

        shell.resize(Nshells);
        Distribute_Mass();  //Define as massas das camadas
        Distribute_Radius();  //Define os raios das camadas
        Distribute_Temperature();
        Distribute_Abundance();
        Initial_RadialVelocity();

        Copy_InitialAbundance();

        reactions = true;
        luminosity = true;
        art_viscosity = true;
        gravity = true;

        eos.id_eos = 2; //ions + photons + pairs

        Update_Interior();
        
        dt_integrator = 0.001;
        
        using PhysConstants::M_Sun;

        std::cerr << "====== White Dwarf Star =======" << std::endl;
        std::cerr << "===============================" << std::endl;
        std::cerr << "Star mass: " << Mass/M_Sun << " M_Sun" << std::endl;
        std::cerr << "Number of shells: " << Nshells << std::endl;
        std::cerr << "===============================" << std::endl;
    }

    Star(){
        std::cerr << "Creating the stellar interior from file...done." << std::endl;

        ReadConfigFile();

        shell.resize(Nshells);

        Read_InitialCondition();
        
        Copy_InitialAbundance();

        Update_Interior();
        
        dt_integrator = 0.001;
    }
 
    ~Star(){}

    void SetMass(const double &m);
    void SetTemperature(const double &T);

    void Distribute_Mass();
    void Distribute_Radius();
    void Distribute_Temperature();
    void Distribute_Abundance();
    void Initial_RadialVelocity();
    void Update_Mass();
    
    void Copy_InitialAbundance();

    void Update_Interior();

    double Kinetic_Energy();
    double Gravitational_Energy();
    double Internal_Energy();
    double Nuclear_Energy();

    void Update_Radius();
    void Update_Velocity();
    void Calculate_Luminosity();

    void dRdt();
    void dRdotdt();
    void ThermalConductionEquation();
    void ThermalConductionEquationDelayed();
    void dQddt();

    void NuclearReactions(double &dt, const double &rtol, const double &tol);
    
    void Single_Step(const double &dt, const double &atol, const double &rtol);
    
    void Fehlberg_Integrator(const double &ti, double &tf, double &dt, const double &atol, const double &rtol);

    void Time_Integration(double &ti, double &dt, const double &rtol, const double &tol, unsigned int id_dynamics = 3);
    
    void HydrostaticEquilibrium(const double &rtol);

    double TimeStep(unsigned int id_dynamics);
    
    
private:

    void T_matrix(std::vector<double>& a, std::vector<double>& b);
    void Q_matrix(std::vector<double>& a, std::vector<double>& b);
    void F_vector(std::vector<double>& f);

    void ReadConfigFile();
    void ProcessLine(std::string *line);
    void Read_InitialCondition();
    void ReadCompositionFile();
    
};

#endif
