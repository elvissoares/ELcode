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

#include "star.hpp"
#include "user_IO.hpp"
#include "util/nr.hpp"


int main (void) {
    
    std::cerr << "The max number of threads is " << omp_get_max_threads() << std::endl;
   
   	time_t begin, end;

    begin = time(0);
    
    using PhysConstants::yr;
    using PhysConstants::M_Sun;
    using PhysConstants::four_pi_o3;
    
    const double tol = 1.e-8, rtol=1.e-6;

     std::vector<double> X(7);
    X[1] = 0.5;
    X[2] = 0.5;
    
    Star star(580,1.16*M_Sun,5.0e7,1.0e8,X);
    star.tau = 0.0;
    star.reactions = true;
    star.luminosity = true;
    star.art_viscosity = true;

   std::cerr << "Reproducing the Case A model from Nomoto 1982...";

   #pragma omp parallel for
   for (std::vector<int>::size_type a = 0; a < star.shell.size(); a++){
       star.shell[a].T = pow(10.0, 7.53 + 0.47*a/565.0);
       if (a > 540) {
           star.shell[a].X[0] = 1.0;
           star.shell[a].X[1] = 0.0;
           star.shell[a].X[2] = 0.0;
       }
       
   }

   #pragma omp parallel for
   for (std::vector<int>::size_type a = 565; a < star.shell.size(); a++){
       star.shell[a].T = pow(10.0, 8.00 - 0.4*(a-565.0)/14.0);
   }
    
   star.shell[565].T = pow(10.0,9.50);

   star.HydrostaticEquilibrium(rtol);

     //Exportanto perfil inicial
    std::string FileName = star.OutputFilePath + "CaseA-profile-N580.dat";
    std::ofstream OutFile(FileName.c_str());
    OutFile << "# M_r   Rho    T " << std::endl;
    
    for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
        OutFile << (j+1)*star.shell[j].M/M_Sun << " " << star.shell[j].Rho << " " << star.shell[j].T << std::endl;
    }
    OutFile.close();
  
    double dt = 1e3*yr, t = 0.;

    OutDynamics myout(star);
    myout.AllVariables(0.);

//     FileName = star.OutputFilePath + "CaseA-shellhelium.dat";
//     OutFile.open(FileName.c_str());
//     OutFile << "# t(yr)   Rho   T " << std::endl;

// //    while ( dt > 0.001 ){
// //          
// //        star.Time_Integration(t,dt,rtol,tol,2);
// //
// //        OutFile << t/yr << " " << star.shell[1130].Rho << " " << star.shell[1130].T  << std::endl;
// //
// //        myout.AllVariables(t/yr);
// //          
// //        std::cerr << "t =" << t/yr << " yr" << std::endl;
// //        std::cerr << "dt =" << dt/yr << std::endl;
// //    }

//     OutFile.close();
    
    std::cerr << "done!" << std::endl;

    //star.HydrostaticEquilibrium(rtol);

    FileName = star.OutputFilePath + "CaseA-N580-onepeak.dat";
    OutFile.open(FileName.c_str());
    OutFile << "# M    R    Rdot     T      X_He     X_C    X_O    X_Ne    X_Mg    X_Si    X_Ni" << std::endl;
    OutFile << "# Energy = " << star.Gravitational_Energy() + star.Internal_Energy() << " ergs" << std::endl;
    
    for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
        OutFile << star.shell[j].M << " " << star.shell[j].Rout << " " << star.shell[j].Rdot_out << " " << star.shell[j].T << " " << star.shell[j].X[0] << " " << star.shell[j].X[1] << " " << star.shell[j].X[2] << " " << star.shell[j].X[3] << " " << star.shell[j].X[4] << " " << star.shell[j].X[5] << " " << star.shell[j].X[6] << std::endl;
    }
	
	end = time(0);

	std::cerr << "Time elapsed: " << difftime(end,begin) << " s"<< std::endl;

	return 0;
}
