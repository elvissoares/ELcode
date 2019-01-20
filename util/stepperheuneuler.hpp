#ifndef STEPPERHEUNEULER_HPP_
#define STEPPERHEUNEULER_HPP_

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
struct StepperHeunEuler : StepperBase {
//Passo de Runge-Kutta quinta ordem de Dormand-Prince com monitoramento de erro de truncamento local para garantir acurácia e ajuste de tamanho de passo.
	typedef D Dtype;		//Disponibiliza o tipo de derivs para odeint.
	std::vector<double> k2,k3,k4,k5,k6;
	std::vector<double> rcont1,rcont2,rcont3,rcont4,rcont5;
	std::vector<double> dydxnew;
	StepperHeunEuler(std::vector<double> &yy, std::vector<double> &dydxx, double &xx,
		const double atoll, const double rtoll, bool dens);
	void step(const double htry,D &derivs);
	void dy(const double h,D &derivs);
	void prepare_dense(const double h,D &derivs);
	double dense_out(const int i, const double x, const double h);
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
StepperHeunEuler<D>::StepperHeunEuler(std::vector<double> &yy,std::vector<double> &dydxx,double &xx, const double atoll,const double rtoll,bool dens) :
	StepperBase(yy,dydxx,xx,atoll,rtoll,dens), k2(n),k3(n),k4(n),k5(n),k6(n),
	rcont1(n),rcont2(n),rcont3(n),rcont4(n),rcont5(n),dydxnew(n) {
	EPS=std::numeric_limits<double>::epsilon();
}

template <class D>
StepperHeunEuler<D>::Controller::Controller() : reject(false), errold(1.0e-4) {}
template <class D>
bool StepperHeunEuler<D>::Controller::success(const double err,double &h) {
	static const double beta=0.0,alpha=0.2-beta*0.75,safe=0.9,minscale=0.2,
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
double StepperHeunEuler<D>::dense_out(const int i,const double x,const double h) {
	double s=(x-xold)/h;
	double s1=1.0-s;
	return rcont1[i]+s*(rcont2[i]+s1*(rcont3[i]+s*(rcont4[i]+s1*rcont5[i])));
}
template <class D>
void StepperHeunEuler<D>::dy(const double h,D &derivs) {
	static const double a21=1/5., b1=0.5, b2=0.5, db1=-0.5, db2=0.5, c2=1.;
    
    std::vector<double> ytemp(y.size());

    // First Step
    for (std::vector<int>::size_type i = 0; i < y.size(); i++) 
        ytemp[i]=y[i]+a21*h*dydx[i];

    derivs(x+c2*h,ytemp,k2);
    
    // Accumulate increments with proper weights, and estimate error as difference between first- and second-order methods.
   
    for (std::vector<int>::size_type i = 0; i < y.size(); i++){
        yout[i] = y[i]+h*(b1*dydx[i]+b2*k2[i]);
        yerr[i] = h*(db1*dydx[i]+db2*k2[i]);
    }

    double xph=x+h;
	derivs(xph,yout,dydxnew);

    ytemp.clear();
}

template <class D>
double StepperHeunEuler<D>::error() {
	double err=0.0,sk;
	for (int i=0;i<n;i++) {
		sk=atol+rtol*MAX(fabs(y[i]),fabs(yout[i]));
		err += SQR(yerr[i]/sk);
	}
	return sqrt(err/n);
}

template <class D>
void StepperHeunEuler<D>::prepare_dense(const double h,D &derivs) {
	std::vector<double> ytemp(n);
	static const double d1=-12715105075.0/11282082432.0,
	d3=87487479700.0/32700410799.0, d4=-10690763975.0/1880347072.0,
	d5=701980252875.0/199316789632.0, d6=-1453857185.0/822651844.0,
	d7=69997945.0/29380423.0;
	for (int i=0;i<n;i++) {
		rcont1[i]=y[i];
		double ydiff=yout[i]-y[i];
		rcont2[i]=ydiff;
		double bspl=h*dydx[i]-ydiff;
		rcont3[i]=bspl;
		rcont4[i]=ydiff-h*dydxnew[i]-bspl;
		rcont5[i]=h*(d1*dydx[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+
			d7*dydxnew[i]);
	}
}

template <class D>
void StepperHeunEuler<D>::step(const double htry,D &derivs) {
	double h=htry;
	for (;;) {
		dy(h,derivs);
		double err= error();

		if (con.success(err,h)) break;
		if (fabs(h) <= fabs(x)*EPS){
			std::cerr << "ERROR: StepperHeunEuler: stepsize underflow!" << std::endl;
			exit(0);
		}
	}
	if (dense)
		prepare_dense(h,derivs);
	dydx=dydxnew;
	y=yout;
	xold=x;
	x += (hdid=h);
	hnext=con.hnext;
}

#endif