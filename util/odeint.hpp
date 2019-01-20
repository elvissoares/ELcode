#ifndef ODEINT_HPP_
#define ODEINT_HPP_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits> 

#include "nr.hpp"
#include "stepper.hpp"
#include "steppersie.hpp"
//#include "stepperbs.hpp"
#include "stepperdopr5.hpp"
#include "stepperfehlb.hpp"

template<class Stepper>
struct Odeint {
	static const int MAXSTP=50000;
	double EPS;
	int nok;
	int nbad;
	int nvar;
	double x1,x2,hmin;
	std::vector<double> y,dydx;
	std::vector<double> &ystart;
	typename Stepper::Dtype &derivs;
	Stepper s;
	int nstp;
	double x,h;
	Odeint(std::vector<double> &ystartt,const double xx1,const double xx2, const double atol,const double rtol,const double h1, const double hminn,typename Stepper::Dtype &derivss);
	void integrate();
};

template<class Stepper>
Odeint<Stepper>::Odeint(std::vector<double> &ystartt, const double xx1, const double xx2, const double atol, const double rtol, const double h1, const double hminn,typename Stepper::Dtype &derivss) : nvar(ystartt.size()), y(nvar),dydx(nvar),ystart(ystartt),x(xx1),nok(0),nbad(0), x1(xx1),x2(xx2),hmin(hminn),derivs(derivss), s(y,dydx,x,atol,rtol) {

	EPS=std::numeric_limits<double>::epsilon();
	h=SIGN(h1,x2-x1);
	for (int i=0;i<nvar;i++) y[i]=ystart[i];
}

template<class Stepper>
void Odeint<Stepper>::integrate() {
	derivs(x,y,dydx);

	for (nstp=0;nstp<MAXSTP;nstp++) {
		if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
			h=x2-x;
		s.step(h,derivs);
		if (s.hdid == h) ++nok; else ++nbad;
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (int i=0;i<nvar;i++) ystart[i]=y[i];
			return;
		}
		if (fabs(s.hnext) <= hmin) {
			std::cerr << "ERROR: Odeint: Step size too small!" << std::endl;
			exit(0);
		}
		h=s.hnext;
	}
	std::cerr << "ERROR: Odeint: Too many steps!" << std::endl;
	exit(0);
}

#endif
