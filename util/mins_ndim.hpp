#ifndef MINS_NDIM_HPP
#define MINS_NDIM_HPP

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <stdexcept>

#include "mins.hpp"

template <class T>
struct F1dim {
	const std::vector<double> &p;
	const std::vector<double> &xi;
	int n;
	T &func;
	std::vector<double> xt;
	F1dim(std::vector<double> &pp, std::vector<double> &xii, T &funcc) : p(pp), xi(xii), n(pp.size()), func(funcc), xt(n) {}
	double operator() (const double x)
	{
		for (int j=0;j<n;j++)
			xt[j]=p[j]+x*xi[j];
		return func(xt);
	}
};
template <class T>
struct Linemethod {
	std::vector<double> p;
	std::vector<double> xi;
	T &func;
	int n;
	Linemethod(T &funcc) : func(funcc) {}
	double linmin()
	{
		double ax,xx,xmin;
		n=p.size();
		F1dim<T> f1dim(p,xi,func);
		ax=0.0;
		xx=1.0;
		Brent brent;
		brent.bracket(ax,xx,f1dim);
		xmin=brent.minimize(f1dim);
		for (int j=0;j<n;j++) {
			xi[j] *= xmin;
			p[j] += xi[j];
		}
		return brent.fmin;
	}
};
template <class T>
struct Df1dim {
	const std::vector<double> &p;
	const std::vector<double> &xi;
	int n;
	T &funcd;
	std::vector<double> xt;
	std::vector<double> dft;
	Df1dim(std::vector<double> &pp, std::vector<double> &xii, T &funcdd) : p(pp), xi(xii), n(pp.size()), funcd(funcdd), xt(n), dft(n) {}
	double operator()(const double x)
	{
		for (int j=0;j<n;j++)
			xt[j]=p[j]+x*xi[j];
		return funcd(xt);
	}
	double df(const double x)
	{
		double df1=0.0;
		funcd.df(xt,dft);
		for (int j=0;j<n;j++)
			df1 += dft[j]*xi[j];
		return df1;
	}
};
template <class T>
struct Dlinemethod {
	std::vector<double> p;
	std::vector<double> xi;
	T &func;
	int n;
	Dlinemethod(T &funcc) : func(funcc) {}
	double linmin()
	{
		double ax,xx,xmin;
		n=p.size();
		Df1dim<T> df1dim(p,xi,func);
		ax=0.0;
		xx=1.0;
		Dbrent dbrent;
		dbrent.bracket(ax,xx,df1dim);
		xmin=dbrent.minimize(df1dim);
		for (int j=0;j<n;j++) {
			xi[j] *= xmin;
			p[j] += xi[j];
		}
		return dbrent.fmin;
	}
};
template <class T>
struct Powell : Linemethod<T> {
	int iter;
	double fret;
	using Linemethod<T>::func;
	using Linemethod<T>::linmin;
	using Linemethod<T>::p;
	using Linemethod<T>::xi;
	const double ftol;
	Powell(T &func, const double ftoll=3.0e-8) : Linemethod<T>(func), ftol(ftoll) {}
	std::vector<double> minimize(std::vector<double> &pp)
	{
		int n=pp.size();
		std::vector< std::vector<double> > ximat(n);
		for (unsigned int i=0;i<n;i++) ximat[i].resize(n);

		for (int i=0;i<n;i++) ximat[i][i]=1.0;
			
		return minimize(pp,ximat);
	}
	std::vector<double> minimize(std::vector<double> &pp, std::vector< std::vector<double> > &ximat)
	{
		const int ITMAX=200;
		const double TINY=1.0e-25;
		double fptt;
		int n=pp.size();
		p=pp;
		std::vector<double> pt(n),ptt(n);
		xi.resize(n);
		fret=func(p);
		for (int j=0;j<n;j++) pt[j]=p[j];
		for (iter=0;;++iter) {
			double fp=fret;
			int ibig=0;
			double del=0.0;
			for (int i=0;i<n;i++) {
				for (int j=0;j<n;j++) xi[j]=ximat[j][i];
				fptt=fret;
				fret=linmin();
				if (fptt-fret > del) {
					del=fptt-fret;
					ibig=i+1;
				}
			}
			if (2.0*(fp-fret) <= ftol*(fabs(fp)+fabs(fret))+TINY) {
				return p;
			}
			if (iter == ITMAX) {
				std::cerr << "ERROR: Powell: exceeding maximum iterations." << std::endl;
				exit(0);
			}
			for (int j=0;j<n;j++) {
				ptt[j]=2.0*p[j]-pt[j];
				xi[j]=p[j]-pt[j];
				pt[j]=p[j];
			}
			fptt=func(ptt);
			if (fptt < fp) {
				double t=2.0*(fp-2.0*fret+fptt)*SQR(fp-fret-del)-del*SQR(fp-fptt);
				if (t < 0.0) {
					fret=linmin();
					for (int j=0;j<n;j++) {
						ximat[j][ibig-1]=ximat[j][n-1];
						ximat[j][n-1]=xi[j];
					}
				}
			}
		}
	}
};
template <class T>
struct Frprmn : Linemethod<T> {
	int iter;
	double fret;
	using Linemethod<T>::func;
	using Linemethod<T>::linmin;
	using Linemethod<T>::p;
	using Linemethod<T>::xi;
	const double ftol;
	Frprmn(T &funcd, const double ftoll=3.0e-8) : Linemethod<T>(funcd), ftol(ftoll) {}
	std::vector<double> minimize(std::vector<double> &pp)
	{
		const int ITMAX=200;
		const double EPS=1.0e-18;
		const double GTOL=1.0e-8;
		double gg,dgg;
		int n=pp.size();
		p=pp;
		std::vector<double> g(n),h(n);
		xi.resize(n);
		double fp=func(p);
		func.df(p,xi);
		for (int j=0;j<n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j];
		}
		for (int its=0;its<ITMAX;its++) {
			iter=its;
			fret=linmin();
			if (2.0*fabs(fret-fp) <= ftol*(fabs(fret)+fabs(fp)+EPS))
				return p;
			fp=fret;
			func.df(p,xi);
			double test=0.0;
			double den=MAX(fabs(fp),1.0);
			for (int j=0;j<n;j++) {
				double temp=fabs(xi[j])*MAX(fabs(p[j]),1.0)/den;
				if (temp > test) test=temp;
			}
			if (test < GTOL) return p;
			dgg=gg=0.0;
			for (int j=0;j<n;j++) {
				gg += g[j]*g[j];
//			  dgg += xi[j]*xi[j];
				dgg += (xi[j]+g[j])*xi[j];
			}
			if (gg == 0.0)
				return p;
			double gam=dgg/gg;
			for (int j=0;j<n;j++) {
				g[j] = -xi[j];
				xi[j]=h[j]=g[j]+gam*h[j];
			}
		}
		throw("Too many iterations in frprmn");
	}
};

#endif
