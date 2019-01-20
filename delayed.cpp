#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <time.h>
#include <math.h>
#include <sstream>
#include <fstream>

#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 0
#define omp_get_thread_num() 0
#endif

#include "star.hpp"
#include "user_IO.hpp"
#include "util/nr.hpp"

template <typename T>
std::string to_string(T const& value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}

int main(int argc, char *argv[]){
    
    std::cerr << "The max number of threads is " << omp_get_max_threads() << std::endl;
   
   	time_t begin, end;

    begin = time(0);
    
    double tau;
    
    if (argc > 1){
        tau = atof(argv[argc-1]);
        std::cout << "The time delay in ms: " << tau << std::endl;
    }
    else {
        std::cout << "Enter the time delay in ms: ";
        std::cin >> tau;
    }
    
    using PhysConstants::yr;
    using PhysConstants::M_Sun;
    using PhysConstants::four_pi_o3;
    
    Star star;

    star.tau = tau/1e3;

    star.OutputFilePath = star.OutputFilePath + "tau" + to_string(tau) + "ms/";

    std::cout << "Analyzing the output directory...";

    std::string DirectoryCommand = "mkdir -p " + star.OutputFilePath;

    system(DirectoryCommand.c_str());

    std::cout << "done!" << std::endl;
  
    const double tol = 1.e-6, rtol=1.e-4;
  
    double dt = 1.e-4, t = 0., tlast=t;
    
    OutDynamics myout(star);
    myout.AllVariables(0.);

    while (t < 2.5){
          
        star.Time_Integration(t,dt,rtol,tol,3);

        if (t - tlast > 0.001) {

            tlast = t;

            myout.AllVariables(t);
          
            std::cerr << "t =" << t << " s" << std::endl;
            std::cerr << "------------------" << std::endl;
        }
    }
    
    std::cerr << "done!" << std::endl;
	
	end = time(0);

	std::cerr << "Time elapsed: " << difftime(end,begin) << " s"<< std::endl;

	return 0;
}
