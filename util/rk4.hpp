#ifndef RK4_HPP_
#define RK4_HPP_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>

#include "nr.hpp"

void rk4(std::vector<double> &y, std::vector<double> &dydx, const double x, const double h, std::vector<double> &yout, void derivs(const double, std::vector<double> &, std::vector<double> &))
{
	int n=y.size();
	std::vector<double> dym(n),dyt(n),yt(n);
	double hh=h*0.5;
	double h6=h/6.0;
	double xh=x+hh;
	for (int i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
	derivs(xh,yt,dyt);
	for (int i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
	derivs(xh,yt,dym);
	for (int i=0;i<n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	derivs(x+h,yt,dyt);
	for (int i=0;i<n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
}

#endif
