#ifndef _CONDUCTIVITY_HPP_
#define _CONDUCTIVITY_HPP_

#include <cmath>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

#include "util/quadrature.hpp"

struct ThermalConductivity {
public:
	double Evaluate_ThermalConductivity(const double &Rhoin, const double &Tin, const std::vector<double> &Yin);
	double Nondegenerateconductivity(const double &Rhoin, const double &Tin, const std::vector<double> &Yin);
	double Degenerateconductivity(const double &logRho, const double &logT, const std::vector<double> &Y);

	ThermalConductivity(){
		Read_Table();
	}

	std::vector<double> RHO, T, Z;
    std::vector< std::vector< std::vector<double> > > KAP;

private:
	void Read_Table();
	void Num_Derivative(std::vector<double> &y, std::vector<double> &y1, std::vector<double> &y2, std::vector<double> &y12, const int &iZ, const int &iR, const int &iT);
};

struct ExponentialIntegral {

	double operator() (const double &t){
		//Integrating with the transformation t=exp(-x)
		return -1.0/log(t);
	}
};

#endif

