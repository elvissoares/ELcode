#ifndef FERMI_HPP_
#define FERMI_HPP_

#include <cmath>
#include <stdio.h>
#include <iostream>

struct FermiDirac {

	double kk, etaa, thetaa;
	double operator() (const double &t);
	double operator() (const double &x, const double &del);
	double F(const double &k, const double &eta, const double &theta);
};

struct FermiDiracDBeta {

	double kk, etaa, thetaa;
	double operator() (const double &t);
	double operator() (const double &x, const double &del);
	double dFdbeta(const double &k, const double &eta, const double &theta);
};

struct FermiDiracDEta {

	double kk, etaa, thetaa;
	double operator() (const double &t);
	double operator() (const double &x, const double &del);
	double dFdeta(const double &k, const double &eta, const double &theta);
};

#endif