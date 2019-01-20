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

#include "ELcode/util/nr.hpp"
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
    
    std::cerr << "================== Accreting mass in a White Dwarf ================== " << std::endl;
    std::cerr << "We are reproducing the Case A from (Nomoto ApJ 1982)" << std::endl;
    
    using PhysConstants::M_Sun;
    using PhysConstants::R_Sun;
    using PhysConstants::G;
    using PhysConstants::four_pi;
    using PhysConstants::four_pi_o3;
    using PhysConstants::yr;
    
    //Reading the pre-existent file with the WD cooled from 10^8 years
    Star star;

    star.OutputFilePath = "output/accretion/";
  
    const double tol = 1.e-6, rtol=1.e-4;
    double dt, t = 0., tlast = 0.; bool hydroeq = true;

    double Mdot = 3.e-8*M_Sun/yr;

    double dtacc = 0.25*star.shell[0].M/Mdot;
    dt = dtacc;

    std::vector<double> Eantes(150);
    std::vector<double> Edepois(150);
    
    std::string FileName;;
    std::ofstream OutFile;

    std::string FileName0 = star.OutputFilePath + "shell-accreted.dat";
    std::ofstream OutFile0(FileName0.c_str());

    std::string FileName2 = star.OutputFilePath + "shell-central.dat";
    std::ofstream OutFile2(FileName2.c_str());

    OutFile0 << 0 << " " << log10(star.shell[53].Rho) << " " << log10(star.shell[53].T) << std::endl;
    OutFile2 << 0 << " " << log10(star.shell[0].Rho) << " " << log10(star.shell[0].T) << std::endl;
    
    double DeltaM = 0.;
    double sumMass;

    while (t < 1e7*yr){
        
        star.shell[star.Nshells-1].M += Mdot*dt;
        DeltaM += Mdot*dt;

        if (star.shell[star.Nshells-1].M > 2.0*star.shell[0].M) {
            
            //Se a massa da camada for maior que um certo limite, adicionamos uma nova camada a estrela

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

            star.shell[star.Nshells-2].Rout = pow(C(star.shell[star.Nshells-2].Rin) + star.shell[star.Nshells-2].M/(four_pi_o3*star.shell[star.Nshells-2].Rho),1/3.);

            star.shell[star.Nshells-1].Rin = star.shell[star.Nshells-2].Rout;
            star.shell[star.Nshells-1].Rout = pow(C(star.shell[star.Nshells-1].Rin) + star.shell[star.Nshells-1].M/(four_pi_o3*star.shell[star.Nshells-2].Rho),1/3.);
        
            star.shell[star.Nshells-1].Rdot_in = star.shell[star.Nshells-2].Rdot_out;
            star.shell[star.Nshells-1].Rdot_out = 0.;
            
            star.Update_Interior();
            
            FileName = star.OutputFilePath + "WDprofile-accreted-t" + to_string(t/yr) + "-Mdot" + to_string(Mdot/(M_Sun/yr) ) + ".dat";
            OutFile.open(FileName.c_str());
            OutFile << "# M(M_sun)    Rho     T" << std::endl;
            
            sumMass = 0.;
            for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
                sumMass += star.shell[j].M;
                OutFile << sumMass/M_Sun << " " << star.shell[j].Rho << " " << star.shell[j].T << std::endl;
            }

            OutFile.close();

            #pragma omp parallel for
            for (std::vector<int>::size_type j = 0; j < star.Nshells; j++)
                Eantes[j] = star.shell[j].Vg + star.shell[j].U;
            
        }
        
        else{
            
            star.shell[star.Nshells-1].Rin = star.shell[star.Nshells-2].Rout;
            star.shell[star.Nshells-1].Rout = pow(C(star.shell[star.Nshells-1].Rin) + star.shell[star.Nshells-1].M/(four_pi_o3*star.shell[star.Nshells-1].Rho),1/3.);
            
            star.Update_Interior();
            
            std::cerr << "star.shell[star.Nshells-1].T: " << star.shell[star.Nshells-1].T << std::endl;

          #pragma omp parallel for
            for (std::vector<int>::size_type j = 0; j < star.Nshells; j++)
                Eantes[j] = star.shell[j].Vg + star.shell[j].U;

            //Continuamos aumentando a massa da camada e contabilizando a mudança das variáveis mecânicas
            star.HydrostaticEquilibrium(rtol);
            
        }
        
        std::cerr << "Delta M = " << DeltaM/M_Sun << " Msun" << std::endl;
        
       #pragma omp parallel for
        for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
            Edepois[j] = star.shell[j].Vg + star.shell[j].U;
            star.shell[j].dUdt = (Eantes[j] - Edepois[j])/dt;
        }

        std::cerr << "dUdt = "<< star.shell[star.Nshells-1].dUdt  << std::endl;
        std::cerr << "L = "<< star.shell[star.Nshells-1].L  << std::endl;

        star.Time_Integration(t,dt,rtol,tol,2);

        dt = dtacc;
        
        std::cerr << "t = " << t/yr << " yr" << std::endl;

    
        if (star.Nshells == 54) OutFile0 << t/yr << " " << log10(star.shell[53].Rho) << " " << log10(star.shell[53].T) << std::endl;
        else OutFile0 << t/yr << " " << log10(star.shell[54].Rho) << " " << log10(star.shell[54].T) << std::endl;
        OutFile2 << t/yr << " " << log10(star.shell[0].Rho) << " " << log10(star.shell[0].T) << std::endl;

    } 

    FileName = star.OutputFilePath + "WD-accreted-DeltaM" + to_string(DeltaM/M_Sun) + "-CaseA.dat";
    OutFile.open(FileName.c_str());
    OutFile << "# M    R    Rdot     T      Y_He     Y_C    Y_O    Y_Ne    Y_Mg    Y_Si    Y_Ni" << std::endl;
    OutFile << "# Energy = " << star.Gravitational_Energy() + star.Internal_Energy() << " ergs" << std::endl;
    
    for (std::vector<int>::size_type j = 0; j < star.Nshells; j++){
        OutFile << star.shell[j].M << " " << star.shell[j].Rout << " " << star.shell[j].Rdot_out << " " << star.shell[j].T << " " << star.shell[j].Y[0] << " " << star.shell[j].Y[1] << " " << star.shell[j].Y[2] << " " << star.shell[j].Y[3] << " " << star.shell[j].Y[4] << " " << star.shell[j].Y[5] << " " << star.shell[j].Y[6] << std::endl;
    }
    
    std::cerr << "done!" << std::endl;
	
	end = time(0);

	std::cerr << "Time elapsed: " << difftime(end,begin) << " s"<< std::endl;

	return 0;
}
