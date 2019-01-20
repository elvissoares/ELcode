#ifndef LUDCMP_HPP_
#define LUDCMP_HPP_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>

#include "nr.hpp"

struct LUdcmp
{
	int n;
	std::vector< std::vector<double> > lu;
	std::vector<int> indx;
	double d;
	LUdcmp(std::vector< std::vector<double> > &a);
	void solve(std::vector<double> &b, std::vector<double> &x);
	void solve(std::vector< std::vector<double> > &b, std::vector< std::vector<double> > &x);
	void inverse(std::vector< std::vector<double> > &ainv);
	double det();
	void mprove(std::vector<double> &b, std::vector<double> &x);
	std::vector< std::vector<double> > &aref;
};
inline LUdcmp::LUdcmp(std::vector< std::vector<double> > &a) : n(a.size()), lu(a), aref(a), indx(n) {
	const double TINY=1.0e-40;
	int i,imax,j,k;
	double big,temp;
	std::vector<double> vv(n);
	d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(lu[i][j])) > big) big=temp;
		if (big == 0.0) std::cerr << "Singular matrix in LUdcmp" << std::endl;
		vv[i]=1.0/big;
	}
	for (k=0;k<n;k++) {
		big=0.0;
		imax=k;
		for (i=k;i<n;i++) {
			temp=vv[i]*fabs(lu[i][k]);
			if (temp > big) {
				big=temp;
				imax=i;
			}
		}
		if (k != imax) {
			for (j=0;j<n;j++) {
				temp=lu[imax][j];
				lu[imax][j]=lu[k][j];
				lu[k][j]=temp;
			}
			d = -d;
			vv[imax]=vv[k];
		}
		indx[k]=imax;
		if (lu[k][k] == 0.0) lu[k][k]=TINY;
		for (i=k+1;i<n;i++) {
			temp=lu[i][k] /= lu[k][k];
			for (j=k+1;j<n;j++)
				lu[i][j] -= temp*lu[k][j];
		}
	}
}
inline void LUdcmp::solve(std::vector<double> &b, std::vector<double> &x)
{
	int i,ii=0,ip,j;
	double sum;
	if (b.size() != n || x.size() != n)
		std::cerr << "LUdcmp::solve bad sizes" << std::endl;
	for (i=0;i<n;i++) x[i] = b[i];
	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=x[ip];
		x[ip]=x[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= lu[i][j]*x[j];
		else if (sum != 0.0)
			ii=i+1;
		x[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=x[i];
		for (j=i+1;j<n;j++) sum -= lu[i][j]*x[j];
		x[i]=sum/lu[i][i];
	}
}

inline void LUdcmp::solve(std::vector< std::vector<double> > &b, std::vector< std::vector<double> > &x)
{
	int i,j,m=b[0].size();
	if (b.size() != n || x.size() != n || b[0].size() != x[0].size())
        std::cerr << "LUdcmp::solve bad sizes" << std::endl;
	std::vector<double> xx(n);
	for (j=0;j<m;j++) {
		for (i=0;i<n;i++) xx[i] = b[i][j];
		solve(xx,xx);
		for (i=0;i<n;i++) x[i][j] = xx[i];
	}
}
inline void LUdcmp::inverse(std::vector< std::vector<double> > &ainv)
{
	ainv.resize(n);
	for (unsigned int i=0;i<n;i++) {
        ainv[i].resize(n);
	}
	solve(ainv,ainv);
}
inline double LUdcmp::det()
{
	double dd = d;
	for (int i=0;i<n;i++) dd *= lu[i][i];
	return dd;
}
inline void LUdcmp::mprove(std::vector<double> &b, std::vector<double> &x)
{
	int i,j;
	std::vector<double> r(n);
	for (i=0;i<n;i++) {
		long double sdp = -b[i];
		for (j=0;j<n;j++)
			sdp += (long double)aref[i][j] * (long double)x[j];
		r[i]=sdp;
	}
	solve(r,r);
	for (i=0;i<n;i++) x[i] -= r[i];
}

#endif
