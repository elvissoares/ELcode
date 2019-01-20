#ifndef ROOTS_HPP
#define ROOTS_HPP

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits> 

#include "nr.hpp"

template <class T>
bool zbrac(T &func, double &x1, double &x2)
{
	const int NTRY=50;
	const double FACTOR=1.6;
	if (x1 == x2) throw("Bad initial range in zbrac");
	double f1=func(x1);
	double f2=func(x2);
	for (int j=0;j<NTRY;j++) {
		if (f1*f2 < 0.0) return true;
		if (fabs(f1) < fabs(f2))
			f1=func(x1 += FACTOR*(x1-x2));
		else
			f2=func(x2 += FACTOR*(x2-x1));
	}
	return false;
}
template <class T>
void zbrak(T &fx, const double x1, const double x2, const int n, std::vector<double> &xb1,
	std::vector<double> &xb2, int &nroot)
{
	int nb=20;
	xb1.resize(nb);
	xb2.resize(nb);
	nroot=0;
	double dx=(x2-x1)/n;
	double x=x1;
	double fp=fx(x1);
	for (int i=0;i<n;i++) {
		double fc=fx(x += dx);
		if (fc*fp <= 0.0) {
			xb1[nroot]=x-dx;
			xb2[nroot++]=x;
			if(nroot == nb) {
				std::vector<double> tempvec1(xb1),tempvec2(xb2);
				xb1.resize(2*nb);
				xb2.resize(2*nb);
				for (int j=0; j<nb; j++) {
					xb1[j]=tempvec1[j];
					xb2[j]=tempvec2[j];
				}
				nb *= 2;
			}
		}
		fp=fc;
	}
}
template <class T>
double rtbis(T &func, const double x1, const double x2, const double xacc) {
	const int JMAX=50;
	double dx,xmid,rtb;
	double f=func(x1);
	double fmid=func(x2);
	if (f*fmid >= 0.0) throw("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (int j=0;j<JMAX;j++) {
		fmid=func(xmid=rtb+(dx *= 0.5));
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	throw("Too many bisections in rtbis");
}
template <class T>
double rtflsp(T &func, const double x1, const double x2, const double xacc) {
	const int MAXIT=30;
	double xl,xh,del;
	double fl=func(x1);
	double fh=func(x2);
	if (fl*fh > 0.0) throw("Root must be bracketed in rtflsp");
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xl=x2;
		xh=x1;
		SWAP(fl,fh);
	}
	double dx=xh-xl;
	for (int j=0;j<MAXIT;j++) {
		double rtf=xl+dx*fl/(fl-fh);
		double f=func(rtf);
		if (f < 0.0) {
			del=xl-rtf;
			xl=rtf;
			fl=f;
		} else {
			del=xh-rtf;
			xh=rtf;
			fh=f;
		}
		dx=xh-xl;
		if (fabs(del) < xacc || f == 0.0) return rtf;
	}
	throw("Maximum number of iterations exceeded in rtflsp");
}
template <class T>
double rtsec(T &func, const double x1, const double x2, const double xacc) {
	const int MAXIT=30;
	double xl,rts;
	double fl=func(x1);
	double f=func(x2);
	if (fabs(fl) < fabs(f)) {
		rts=x1;
		xl=x2;
		SWAP(fl,f);
	} else {
		xl=x1;
		rts=x2;
	}
	for (int j=0;j<MAXIT;j++) {
		double dx=(xl-rts)*f/(f-fl);
		xl=rts;
		fl=f;
		rts += dx;
		f=func(rts);
		if (fabs(dx) < xacc || f == 0.0) return rts;
	}
	throw("Maximum number of iterations exceeded in rtsec");
}
template <class T>
double zriddr(T &func, const double x1, const double x2, const double xacc) {
	const int MAXIT=60;
	double fl=func(x1);
	double fh=func(x2);
	if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		double xl=x1;
		double xh=x2;
		double ans=-9.99e99;
		for (int j=0;j<MAXIT;j++) {
			double xm=0.5*(xl+xh);
			double fm=func(xm);
			double s=sqrt(fm*fm-fl*fh);
			if (s == 0.0) return ans;
			double xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
			if (fabs(xnew-ans) <= xacc) return ans;
			ans=xnew;
			double fnew=func(ans);
			if (fnew == 0.0) return ans;
			if (SIGN(fm,fnew) != fm) {
				xl=xm;
				fl=fm;
				xh=ans;
				fh=fnew;
			} else if (SIGN(fl,fnew) != fl) {
				xh=ans;
				fh=fnew;
			} else if (SIGN(fh,fnew) != fh) {
				xl=ans;
				fl=fnew;
			} else throw("never get here.");
			if (fabs(xh-xl) <= xacc) return ans;
		}
		throw("zriddr exceed maximum iterations");
	}
	else {
		if (fl == 0.0) return x1;
		if (fh == 0.0) return x2;
		throw("root must be bracketed in zriddr.");
	}
}
template <class T>
double zbrent(T &func, const double x1, const double x2, const double tol)
{
	const int ITMAX=100;
	const double EPS=std::numeric_limits<double>::epsilon();
	double a=x1,b=x2,c=x2,d,e,fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
		std::cerr << "Root must be bracketed in zbrent" << std::endl;
		std::cerr << "fa = " << fa << std::endl;
		std::cerr << "fb = " << fb << std::endl;
		exit(0);
	}
	fc=fb;
	for (int iter=0;iter<ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			double min1=3.0*xm*q-fabs(tol1*q);
			double min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
			fb=func(b);
	}
	throw("Maximum number of iterations exceeded in zbrent");
}
template <class T>
double rtnewt(T &funcd, const double x1, const double x2, const double xacc) {
	const int JMAX=20;
	double rtn=0.5*(x1+x2);
	for (int j=0;j<JMAX;j++) {
		double f=funcd(rtn);
		double df=funcd.df(rtn);
		double dx=f/df;
		rtn -= dx;
		if ((x1-rtn)*(rtn-x2) < 0.0)
			throw("Jumped out of brackets in rtnewt");
		if (fabs(dx) < xacc) return rtn;
	}
	throw("Maximum number of iterations exceeded in rtnewt");
}
template <class T>
double rtsafe(T &funcd, const double x1, const double x2, const double xacc) {
	const int MAXIT=100;
	double xh,xl;
	double fl=funcd(x1);
	double fh=funcd(x2);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		throw("Root must be bracketed in rtsafe");
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	double rts=0.5*(x1+x2);
	double dxold=fabs(x2-x1);
	double dx=dxold;
	double f=funcd(rts);
	double df=funcd.df(rts);
	for (int j=0;j<MAXIT;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			double temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		f=funcd(rts);
		df=funcd.df(rts);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
	throw("Maximum number of iterations exceeded in rtsafe");
}

#endif