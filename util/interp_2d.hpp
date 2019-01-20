#ifndef INTERP2D_HPP_
#define INTERP2D_HPP_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits> 

#include "nr.hpp"
#include "interp_1d.hpp"
#include "interp_linear.hpp"

struct Bilin_interp {
	int m,n;
	const std::vector< std::vector<double> > &y;
	Linear_interp x1terp, x2terp;

	Bilin_interp(std::vector<double> &x1v, std::vector<double> &x2v, std::vector< std::vector<double> > &ym) : m(x1v.size()), n(x2v.size()), y(ym), x1terp(x1v,x1v), x2terp(x2v,x2v) {}

	double interp(double x1p, double x2p) {
		int i,j;
		double yy, t, u;
		i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
		j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);
		t = (x1p-x1terp.xx[i])/(x1terp.xx[i+1]-x1terp.xx[i]);
		u = (x2p-x2terp.xx[j])/(x2terp.xx[j+1]-x2terp.xx[j]);
		yy = (1.-t)*(1.-u)*y[i][j] + t*(1.-u)*y[i+1][j] + (1.-t)*u*y[i][j+1] + t*u*y[i+1][j+1];
		return yy;
	}
};
struct Poly2D_interp {
	int m,n,mm,nn;
	const std::vector< std::vector<double> > &y;
	std::vector<double> yv;
	Poly_interp x1terp, x2terp;

	Poly2D_interp(std::vector<double> &x1v, std::vector<double> &x2v, std::vector< std::vector<double> > &ym,
		int mp, int np) : m(x1v.size()), n(x2v.size()),
		mm(mp), nn(np), y(ym), yv(m),
		x1terp(x1v,yv,mm), x2terp(x2v,x2v,nn) {}

	double interp(double x1p, double x2p) {
		int i,j,k;
		i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
		j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);
		for (k=i;k<i+mm;k++) {
			x2terp.yy = &y[k][0];
			yv[k] = x2terp.rawinterp(j,x2p);
		}
		return x1terp.rawinterp(i,x1p);
	}
};
struct Spline2D_interp {
	int m,n;
	std::vector< std::vector<double> > y;
	std::vector<double> x1;
	std::vector<double> yv;
	std::vector<Spline_interp*> srp;

	Spline2D_interp(std::vector<double> &x1v, std::vector<double> &x2v, std::vector< std::vector<double> > &ym)
		: m(x1v.size()), n(x2v.size()), y(ym), yv(m), x1(x1v), srp(m) {
		for (int i=0;i<m;i++) srp[i] = new Spline_interp(x2v,&y[i][0]);
	}

	~Spline2D_interp(){
		srp.clear();
	}

	double interp(double x1p, double x2p) {
		for (int i=0;i<m;i++) yv[i] = (*srp[i]).interp(x2p);
		Spline_interp scol(x1,yv);
		return scol.interp(x1p);
	}
};
void bcucof(std::vector<double> &y, std::vector<double> &y1, std::vector<double> &y2, std::vector<double> &y12, const double d1, const double d2, std::vector< std::vector<double> > &cmat) {
	static int wt_d[16*16]=
		{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
		-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
		2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
		0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
		0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
		0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
		-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
		9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
		-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
		2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
		-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
		4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1};
	int l,k,j,i;
	double xx,d1d2=d1*d2;
	std::vector<double> cl(16),x(16);

	for (i=0;i<4;i++) {
		x[i]=y[i];
		x[i+4]=y1[i]*d1;
		x[i+8]=y2[i]*d2;
		x[i+12]=y12[i]*d1d2;
	}
	for (i=0;i<16;i++) {
		xx=0.0;
		for (k=0;k<16;k++) xx += wt_d[i*16+k]*x[k];
		cl[i]=xx;
	}
	l=0;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++) cmat[i][j]=cl[l++];

	cl.clear(); x.clear();
}
void bcuint(std::vector<double> &y, std::vector<double> &y1, std::vector<double> &y2, std::vector<double> &y12,
	const double x1l, const double x1u, const double x2l, const double x2u,
	const double x1, const double x2, double &ansy, double &ansy1, double &ansy2) {
	int i;
	double t,u,d1=x1u-x1l,d2=x2u-x2l;

	std::vector< std::vector<double> > cmat(4);
	for (i=0;i<4;i++) cmat[i].resize(4);

	bcucof(y,y1,y2,y12,d1,d2,cmat);
	if (x1u == x1l || x2u == x2l)
		throw("Bad input in routine bcuint");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	ansy=ansy2=ansy1=0.0;
	for (i=3;i>=0;i--) {
		ansy=t*ansy+((cmat[i][3]*u+cmat[i][2])*u+cmat[i][1])*u+cmat[i][0];
		ansy2=t*ansy2+(3.0*cmat[i][3]*u+2.0*cmat[i][2])*u+cmat[i][1];
		ansy1=u*ansy1+(3.0*cmat[3][i]*t+2.0*cmat[2][i])*t+cmat[1][i];
	}
	ansy1 /= d1;
	ansy2 /= d2;
}

#endif
