#ifndef STEPPER_HPP_
#define STEPPER_HPP_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>

#include "nr.hpp"

struct StepperBase {
	double &x;
	double xold;
	std::vector<double> &y, &dydx;
	double atol,rtol;
	double hdid;
	double hnext;
	double EPS;
	int n,neqn;
	std::vector<double> yout, yerr;
	StepperBase(std::vector<double> &yy, std::vector<double> &dydxx, double &xx, const double atoll, const double rtoll) : x(xx),y(yy),dydx(dydxx),atol(atoll), rtol(rtoll),n(y.size()),neqn(n),yout(n),yerr(n) {}
};

#endif
