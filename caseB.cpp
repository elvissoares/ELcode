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
   
   	time_t begin, end;

    begin = time(0);
    
    using PhysConstants::yr;
    using PhysConstants::M_Sun;
    using PhysConstants::four_pi_o3;
    
    const double tol = 1.e-8, rtol=1.e-3;

    std::vector<double> X(7);
    X[1] = 0.5;
    X[2] = 0.5;
    
    std::cerr << "Reproducing the Case B model from Nomoto 1982" << std::endl;
    
    Star star(131,1.31*M_Sun,3.0e8,1.0e7,X);

   //#pragma omp parallel for
   for (std::vector<int>::size_type a = 0; a < star.shell.size(); a++){
       star.shell[a].T = pow(10.0, 7.80 - 0.03*a/109.0);
       if (a > 108) {
           star.shell[a].X[0] = 1.0;
           star.shell[a].X[1] = 0.0;
           star.shell[a].X[2] = 0.0;
       }
   }

    star.shell[130].T = pow(10.0,7.55);
   
    star.shell[109].T = pow(10.0,9.5);
    
    std::cerr << "Obtaining the hydrostatic equilibrium..." << std::endl;
    
    star.Update_Interior();
    star.HydrostaticEquilibrium(rtol);
    
    std::cerr << "done!" << std::endl;
   

    //Exportanto perfil inicial
    std::string FileName = star.OutputFilePath + "CaseB-profile-N131.dat";
    std::ofstream OutFile(FileName.c_str());
    OutFile << "# M_r   Rho    T " << std::endl;
    
    for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
        OutFile << (j+1)*star.shell[j].M/M_Sun << " " << star.shell[j].Rho << " " << star.shell[j].T << std::endl;
    }
    OutFile.close();
    

    FileName = star.OutputFilePath + "CaseB-N131-onepeak.dat";
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
