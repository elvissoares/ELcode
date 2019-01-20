#ifndef STEPPERSIE_HPP_
#define STEPPERSIE_HPP_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>

#include "nr.hpp"
#include "stepper.hpp"
#include "ludcmp.hpp"

template <class D>
struct StepperSie : StepperBase {
	typedef D Dtype;
	static const int KMAXX=12,IMAXX=KMAXX+1;
	int k_targ;
	std::vector<int> nseq;
	std::vector<double> cost;
	std::vector< std::vector<double> > table;
	std::vector< std::vector<double> > dfdy;
	std::vector<double> dfdx;
	double jac_redo;
	bool calcjac;
	double theta;
	std::vector< std::vector<double> > a;
	int kright;
	std::vector< std::vector<double> > coeff;
	std::vector< std::vector<double> > fsave;
	std::vector<double> dens;
	std::vector<double> factrl;
	StepperSie(std::vector<double> &yy, std::vector<double> &dydxx, double &xx, const double atol,
		const double rtol);
	void step(const double htry,D &derivs);
	bool dy(std::vector<double> &y, const double htot, const int k, std::vector<double> &yend,
		int &ipt,std::vector<double> &scale,D &derivs);
	void polyextr(const int k, std::vector< std::vector<double> > &table, std::vector<double> &last);
};
template <class D>
StepperSie<D>::StepperSie(std::vector<double> &yy, std::vector<double> &dydxx, double &xx,
	const double atoll,const double rtoll)
	: StepperBase(yy,dydxx,xx,atoll,rtoll),nseq(IMAXX),cost(IMAXX),dfdx(n),calcjac(false),dens((IMAXX+2)*n),factrl(IMAXX) {
     
    table.resize(KMAXX);
    for (int k=0;k<KMAXX;k++) table[k].resize(n);
    
    coeff.resize(IMAXX);
    for (int k=0;k<IMAXX;k++) coeff[k].resize(IMAXX);
            
    a.resize(n);
    dfdy.resize(n);
    for (int k=0;k<n;k++){
        a[k].resize(n);
        dfdy[k].resize(n);
    }
        
    fsave.resize((IMAXX-1)*(IMAXX+1)/2+2);
    for (int k=0;k<(IMAXX-1)*(IMAXX+1)/2+2;k++){
        fsave[k].resize(n);
    }
    
	static const double costfunc=1.0,costjac=5.0,costlu=1.0,costsolve=1.0;
    EPS=std::numeric_limits<double>::epsilon();
	jac_redo=MIN(1.0e-4,rtol);
	theta=2.0*jac_redo;
	nseq[0]=2;
	nseq[1]=3;
	for (int i=2;i<IMAXX;i++)
		nseq[i]=2*nseq[i-2];
	cost[0]=costjac+costlu+nseq[0]*(costfunc+costsolve);
	for (int k=0;k<KMAXX;k++)
		cost[k+1]=cost[k]+(nseq[k+1]-1)*(costfunc+costsolve)+costlu;
	hnext=-1.0e99;
	double logfact=-log10(rtol+atol)*0.6+0.5;
	k_targ=MAX(1,MIN(KMAXX-1,int(logfact)));
	for (int k=0; k<IMAXX; k++) {
		for (int l=0; l<k; l++) {
			double ratio=double(nseq[k])/nseq[l];
			coeff[k][l]=1.0/(ratio-1.0);
		}
	}
	factrl[0]=1.0;
	for (int k=0; k<IMAXX-1; k++)
		factrl[k+1]=(k+1)*factrl[k];
}
template <class D>
void StepperSie<D>::step(const double htry,D &derivs) {
	const double STEPFAC1=0.6,STEPFAC2=0.93,STEPFAC3=0.1,STEPFAC4=4.0,
		STEPFAC5=0.5,KFAC1=0.7,KFAC2=0.9;
	static bool first_step=true,last_step=false;
	static bool forward,reject=false,prev_reject=false;
	static double errold;
	int i,k;
	double fac,h,hnew,err;
	bool firstk;
	std::vector<double> hopt(IMAXX),work(IMAXX);
	std::vector<double> ysav(n),yseq(n);
	std::vector<double> ymid(n),scale(n);
	work[0]=1.e30;
	h=htry;
	forward = h>0 ? true : false;
	for (i=0;i<n;i++) ysav[i]=y[i];
	if (h != hnext && !first_step) {
		last_step=true;
	}
	if (reject) {
		prev_reject=true;
		last_step=false;
		theta=2.0*jac_redo;
	}
	for (i=0;i<n;i++)
		scale[i]=atol+rtol*fabs(y[i]);
	reject=false;
	firstk=true;
	hnew=fabs(h);
	compute_jac:
	if (theta > jac_redo && !calcjac) {
		derivs.jacobian(x,y,dfdx,dfdy);
		calcjac=true;
	}
	while (firstk || reject) {
		h = forward ? hnew : -hnew;
		firstk=false;
		reject=false;
		if (fabs(h) <= fabs(x)*EPS)
           std::cerr << "step size h = " << h << " underflow in StepperSie" << std::endl;
		int ipt=-1;
		for (k=0; k<=k_targ+1;k++) {
			bool success=dy(ysav,h,k,yseq,ipt,scale,derivs);
			if (!success) {
				reject=true;
				hnew=fabs(h)*STEPFAC5;
				break;
			}
			if (k == 0)
				 y=yseq;
			else
				for (i=0;i<n;i++)
					table[k-1][i]=yseq[i];
			if (k != 0) {
				polyextr(k,table,y);
				err=0.0;
				for (i=0;i<n;i++) {
					scale[i]=atol+rtol*fabs(ysav[i]);
					err+=SQR((y[i]-table[0][i])/scale[i]);
				}
				err=sqrt(err/n);
				if (err > 1.0/EPS || (k > 1 && err >= errold)) {
					reject=true;
					hnew=fabs(h)*STEPFAC5;
					break;
				}
				errold=MAX(4.0*err,1.0);
				double expo=1.0/(k+1);
				double facmin=pow(STEPFAC3,expo);
				if (err == 0.0)
					fac=1.0/facmin;
				else {
					fac=STEPFAC2/pow(err/STEPFAC1,expo);
					fac=MAX(facmin/STEPFAC4,MIN(1.0/facmin,fac));
				}
				hopt[k]=fabs(h*fac);
				work[k]=cost[k]/hopt[k];
				if ((first_step || last_step) && err <= 1.0)
					break;
				if (k == k_targ-1 && !prev_reject && !first_step && !last_step) {
					if (err <= 1.0)
						break;
					else if (err>nseq[k_targ]*nseq[k_targ+1]*4.0) {
						reject=true;
						k_targ=k;
						if (k_targ>1 && work[k-1]<KFAC1*work[k])
							k_targ--;
						hnew=hopt[k_targ];
						break;
					}
				}
				if (k == k_targ) {
					if (err <= 1.0)
						break;
					else if (err>nseq[k+1]*2.0) {
						reject=true;
						if (k_targ>1 && work[k-1]<KFAC1*work[k])
							k_targ--;
						hnew=hopt[k_targ];
						break;
					}
				}
				if (k == k_targ+1) {
					if (err > 1.0) {
						reject=true;
						if (k_targ>1 && work[k_targ-1]<KFAC1*work[k_targ])
							k_targ--;
						hnew=hopt[k_targ];
					}
					break;
				}
			}
		}
		if (reject) {
			prev_reject=true;
			if (!calcjac) {
				theta=2.0*jac_redo;
				goto compute_jac;
			}
		}
	}
	calcjac=false;

	xold=x;
	x+=h;
	hdid=h;
	first_step=false;
	int kopt;
	if (k == 1)
		kopt=2;
	else if (k <= k_targ) {
		kopt=k;
		if (work[k-1] < KFAC1*work[k])
			kopt=k-1;
		else if (work[k] < KFAC2*work[k-1])
			kopt=MIN(k+1,KMAXX-1);
	} else {
		kopt=k-1;
		if (k > 2 && work[k-2] < KFAC1*work[k-1])
			kopt=k-2;
		if (work[k] < KFAC2*work[kopt])
			kopt=MIN(k,KMAXX-1);
	}
	if (prev_reject) {
		k_targ=MIN(kopt,k);
		hnew=MIN(fabs(h),hopt[k_targ]);
		prev_reject=false;
	}
	else {
		if (kopt <= k)
			hnew=hopt[kopt];
		else {
			if (k<k_targ && work[k]<KFAC2*work[k-1])
				hnew=hopt[k]*cost[kopt+1]/cost[k];
			else
				hnew=hopt[k]*cost[kopt]/cost[k];
		}
		k_targ=kopt;
	}
	if (forward)
		hnext=hnew;
	else
		hnext=-hnew;	
}
template <class D>
bool StepperSie<D>::dy(std::vector<double> &y,const double htot,const int k,std::vector<double> &yend,
	int &ipt,std::vector<double> &scale,D &derivs) {
	std::vector<double> del(n),ytemp(n),dytemp(n);
	int nstep=nseq[k];
	double h=htot/nstep;
	for (int i=0;i<n;i++) {
		for (int j=0;j<n;j++) a[i][j] = -dfdy[i][j];
		a[i][i] += 1.0/h;
	}
	LUdcmp alu(a);
	double xnew=x+h;
	derivs(xnew,y,del);
	for (int i=0;i<n;i++)
		ytemp[i]=y[i];
	alu.solve(del,del);
	if ( nstep==k+1) {
		ipt++;
		for (int i=0;i<n;i++)
			fsave[ipt][i]=del[i];
	}
	for (int nn=1;nn<nstep;nn++) {
		for (int i=0;i<n;i++)
			ytemp[i] += del[i];
		xnew += h;
		derivs(xnew,ytemp,yend);
		if (nn ==1 && k<=1) {
			double del1=0.0;
			for (int i=0;i<n;i++)
				del1 += SQR(del[i]/scale[i]);
			del1=sqrt(del1);
			derivs(x+h,ytemp,dytemp);
			for (int i=0;i<n;i++)
				del[i]=dytemp[i]-del[i]/h;
			alu.solve(del,del);
			double del2=0.0;
			for (int i=0;i<n;i++)
				del2 += SQR(del[i]/scale[i]);
			del2=sqrt(del2);
			theta=del2/MAX(1.0,del1);
			if (theta > 1.0)
				return false;
		}
		alu.solve(yend,del);
		if ( nn >= nstep-k-1) {
			ipt++;
			for (int i=0;i<n;i++)
				fsave[ipt][i]=del[i];
		}
	}
	for (int i=0;i<n;i++)
		yend[i]=ytemp[i]+del[i];
	return true;
}
template <class D>
void StepperSie<D>::polyextr(const int k,std::vector< std::vector<double> > &table,std::vector<double> &last) {
	int l=last.size();
	for (int j=k-1; j>0; j--)
		for (int i=0; i<l; i++)
			table[j-1][i]=table[j][i]+coeff[k][j]*(table[j][i]-table[j-1][i]);
	for (int i=0; i<l; i++)
		last[i]=table[0][i]+coeff[k][0]*(table[0][i]-last[i]);
}

#endif
