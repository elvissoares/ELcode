#ifndef STEPPERFEHLB_HPP_
#define STEPPERFEHLB_HPP_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits> 

#include "nr.hpp"
#include "stepper.hpp"

template <class D>
struct StepperFehlb : StepperBase {
//Passo de Runge-Kutta quinta ordem de Dormand-Prince com monitoramento de erro de truncamento local para garantir acurácia e ajuste de tamanho de passo.
	typedef D Dtype;		//Disponibiliza o tipo de derivs para odeint.
	std::vector<double> k1,k2;
	std::vector<double> dydxnew;
	StepperFehlb(std::vector<double> &yy, std::vector<double> &dydxx, double &xx,
		const double atoll, const double rtoll);
	void step(const double htry,D &derivs);
	void dy(const double h,D &derivs);
	double error();
	struct Controller {
		double hnext,errold;
		bool reject;
		Controller();
		bool success(const double err, double &h);
	};
	Controller con;
};

template <class D>
StepperFehlb<D>::StepperFehlb(std::vector<double> &yy,std::vector<double> &dydxx,double &xx, const double atoll,const double rtoll) :
	StepperBase(yy,dydxx,xx,atoll,rtoll), k1(n),k2(n),dydxnew(n) {
	EPS=std::numeric_limits<double>::epsilon();
}

template <class D>
StepperFehlb<D>::Controller::Controller() : reject(false), errold(1.0e-4) {}
template <class D>
bool StepperFehlb<D>::Controller::success(const double err,double &h) {
	static const double beta=0.2,alpha=0.5-beta*0.75,safe=0.8,minscale=0.2,
		maxscale=10.0;
	double scale;
	if (err <= 1.0) {
		if (err == 0.0)
			scale=maxscale;
		else {
			scale=safe*pow(err,-alpha)*pow(errold,beta);
			if (scale<minscale) scale=minscale;
			if (scale>maxscale) scale=maxscale;
		}
		if (reject)
			hnext=h*MIN(scale,1.0);
		else
			hnext=h*scale;
		errold=MAX(err,1.0e-4);
		reject=false;
		return true;
	} else {
		scale=MAX(safe*pow(err,-alpha),minscale);
		h *= scale;
		reject=true;
		return false;
	}
}
template <class D>
void StepperFehlb<D>::dy(const double h,D &derivs) {

	std::vector<double> ytemp(n);
	int i;
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+0.5*h*dydx[i];
	derivs(x+0.5*h,ytemp,k1);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+(1./256.)*h*(dydx[i]+255*k1[i]);
	

	double xph=x+h;
	derivs(xph,ytemp,dydxnew);

	for (i=0;i<n;i++){
		yout[i]=y[i]+(1./512.)*h*(dydx[i]+510*k1[i]+dydxnew[i]);
		yerr[i]=yout[i]-ytemp[i];
	}
		
}

template <class D>
double StepperFehlb<D>::error() {
	double err=0.0,sk;
	for (int i=0;i<n;i++) {
		sk=atol+rtol*MAX(fabs(y[i]),fabs(yout[i]));
		err += SQR(yerr[i]/sk);
	}
	return sqrt(err/n);
}

template <class D>
void StepperFehlb<D>::step(const double htry,D &derivs) {
	double h=htry;
	for (;;) {
		dy(h,derivs);
		double err=error();
		if (con.success(err,h)) break;
		if (fabs(h) <= fabs(x)*EPS){
			std::cerr << "ERROR: StepperFehlb: stepsize underflow!" << std::endl;
			exit(0);
		}
	}
	dydx=dydxnew;
	y=yout;
	xold=x;
	x += (hdid=h);
	hnext=con.hnext;
}

#endif
