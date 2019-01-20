#ifndef INTERPLIN_HPP_
#define INTERPLIN_HPP_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits> 

#include "nr.hpp"
#include "interp_1d.hpp"

struct Linear_interp : Base_interp
{
	Linear_interp(std::vector<double> &xv, std::vector<double> &yv)
		: Base_interp(xv,&yv[0],2)  {}
	double rawinterp(int j, double x) {
		if (xx[j]==xx[j+1]) return yy[j];
		else return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
	}
};

#endif
