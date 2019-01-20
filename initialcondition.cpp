#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <time.h>
#include <math.h>
#include <sstream>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 0
#define omp_get_thread_num() 0
#endif

#include "ELcode/star.hpp"
#include "ELcode/user_IO.hpp"

template <typename T>
std::string to_string(T const& value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}

int main (void) {
    
    std::cerr << "The max number of threads is " << omp_get_max_threads() << std::endl;
   
   	time_t begin, end;

    begin = time(0);
    
    using PhysConstants::M_Sun;
    using PhysConstants::R_Sun;
    using PhysConstants::G;
    using PhysConstants::four_pi;
    using PhysConstants::four_pi_o3;
    using PhysConstants::yr;
    
    const double tol = 1.e-6, rtol=1.e-4;

    std::vector<double> X(7);
    X[1] = 0.5;
    X[2] = 0.5;
    
    double Mass;
    int Nshells;
    
    std::cout << "The initialcondition file creates the WD star" << std::endl;
    std::cout << "Enter the Mass (solar masses) of the WD: ";
    std::cin >> Mass;
    
    std::cout << "Enter the Number of shells to be used: ";
    std::cin >> Nshells;
    
    Star star(Nshells,Mass*M_Sun,5.0e7,1.0e9,X);

    std::cerr << "T = "<< star.shell[star.Nshells-1].T  << std::endl;

    star.Update_Interior();
    
    star.HydrostaticEquilibrium(rtol);
    
    star.OutputFilePath = "input/";

    double tnext, dt, t = 0., tlast = 0.;
    
    dt = 10*yr;
    
    //Cooling down the WD during 10^8 years
    
    std::cerr << "----Cooling down the WD during 10^8 years" << std::endl;
    
    star.reactions = false;
    
    std::stringstream Aux;
    Aux << "-M" << star.Mass/M_Sun << "-N" << star.shell.size();
    
    std::string FileName;;
    std::ofstream OutFile;
    
    FileName = star.OutputFilePath + "WD-centralinformation" + Aux.str() + ".dat";
    OutFile.open(FileName.c_str());
    OutFile << "# t(yr)   Rho(g/cm3)   T(K) " << std::endl;
    
    
    while (t < 1e8*yr){
        
        star.Time_Integration(t,dt,rtol,tol,2);
        
        if (t-tlast > 1e7*yr) {

             std::cerr << "t = " << t/yr << " yr" << std::endl;

            OutFile << t/yr << " " << star.shell[0].Rho << " " << star.shell[0].T << std::endl;
            
            tlast=t;
        }
    }
    OutFile.close();
    
    star.HydrostaticEquilibrium(rtol);
    
    //Exportanto perfil inicial
    FileName = star.OutputFilePath + "WD-profile" + Aux.str() + ".dat";
    OutFile.open(FileName.c_str());
    OutFile << "# Rho  T " << std::endl;
    
    for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
        OutFile << star.shell[j].Rho << " " << star.shell[j].T << std::endl;
    }
    OutFile.close();
    
    FileName = star.OutputFilePath + "WD" + Aux.str() + ".dat";
    OutFile.open(FileName.c_str());
    OutFile << "# M    R    Rdot     T      Y_He     Y_C    Y_O    Y_Ne    Y_Mg    Y_Si    Y_Ni" << std::endl;
    OutFile << "# Energy = " << star.Gravitational_Energy() + star.Internal_Energy() << " ergs" << std::endl;
    
    for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
        OutFile << star.shell[j].M << " " << star.shell[j].Rout << " " << star.shell[j].Rdot_out << " " << star.shell[j].T << " " << star.shell[j].Y[0] << " " << star.shell[j].Y[1] << " " << star.shell[j].Y[2] << " " << star.shell[j].Y[3] << " " << star.shell[j].Y[4] << " " << star.shell[j].Y[5] << " " << star.shell[j].Y[6] << std::endl;
    }
	
	end = time(0);

	std::cerr << "Time elapsed: " << difftime(end,begin) << " s"<< std::endl;

	return 0;
}
