#ifndef STEPPERBS_HPP_
#define STEPPERBS_HPP_

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
struct StepperBS : StepperBase {
	typedef D Dtype;
	static const int KMAXX=8,IMAXX=KMAXX+1;
	int k_targ;
	std::vector<int> nseq;
	std::vector<int> cost;
	std::vector< std::vector<double> >table;
	std::vector<double> dydxnew;
	int mu;
	std::vector< std::vector<double> >coeff;
	std::vector<double> errfac;
	std::vector< std::vector<double> >ysave;
	std::vector< std::vector<double> >fsave;
	std::vector<int> ipoint;
	StepperBS(std::vector<double> &yy, std::vector<double> &dydxx, double &xx, const double atol,
		const double rtol);
	void step(const double htry,D &derivs);
	virtual void dy(const std::vector<double> &y, const double htot, const int k, std::vector<double> &yend,
		int &ipt, D &derivs);
	void polyextr(const int k, std::vector< std::vector<double> > &table, std::vector<double> &last);
};
template <class D>
StepperBS<D>::StepperBS(std::vector<double> &yy,std::vector<double> &dydxx,double &xx,
	const double atoll,const double rtoll) :
	StepperBase(yy,dydxx,xx,atoll,rtoll),nseq(IMAXX),cost(IMAXX),
	dydxnew(n),errfac(2*IMAXX+2),ipoint(IMAXX+1) {
        EPS=std::numeric_limits<double>::epsilon();
        
        table.resize(KMAXX);
        ysave.resize(KMAXX);
        for (int k=0;k<KMAXX;k++) {
            table[k].resize(n);
            ysave[k].resize(n);
        }
        
        coeff.resize(IMAXX);
        for (int k=0;k<IMAXX;k++) coeff[k].resize(IMAXX);
        
        fsave.resize((IMAXX-1)*(IMAXX+1)/2+2);
        for (int k=0;k<(IMAXX-1)*(IMAXX+1)/2+2;k++){
            fsave[k].resize(n);
        }
        
    for (int i=0;i<IMAXX;i++)
        nseq[i]=2*(i+1);
	cost[0]=nseq[0]+1;
	for (int k=0;k<KMAXX;k++) cost[k+1]=cost[k]+nseq[k+1];
	hnext=-1.0e99;
	double logfact=-log10(MAX(1.0e-12,rtol))*0.6+0.5;
	k_targ=MAX(1,MIN(KMAXX-1,int(logfact)));
	for (int k = 0; k<IMAXX; k++) {
		for (int l=0; l<k; l++) {
			double ratio=double(nseq[k])/nseq[l];
			coeff[k][l]=1.0/(ratio*ratio-1.0);
		}
	}
	for (int i=0; i<2*IMAXX+1; i++) {
		int ip5=i+5;
		errfac[i]=1.0/(ip5*ip5);
		double e = 0.5*sqrt(double(i+1)/ip5);
		for (int j=0; j<=i; j++) {
			errfac[i] *= e/(j+1);
		}
	}
	ipoint[0]=0;
	for (int i=1; i<=IMAXX; i++) {
		int njadd=4*i-2;
		if (nseq[i-1] > njadd) njadd++;
		ipoint[i]=ipoint[i-1]+njadd;
	}
}
template <class D>
void StepperBS<D>::dy(const std::vector<double> &y,const double htot,const int k,std::vector<double> &yend,
	int &ipt,D &derivs) {
	std::vector<double> ym(n),yn(n);
	int nstep=nseq[k];
	double h=htot/nstep;
	for (int i=0;i<n;i++) {
		ym[i]=y[i];
		yn[i]=y[i]+h*dydx[i];
	}
	double xnew=x+h;
	derivs(xnew,yn,yend);
	double h2=2.0*h;
	for (int nn=1;nn<nstep;nn++) {
		if (nn == nstep/2) {
				for (int i=0;i<n;i++)
					ysave[k][i]=yn[i];
		}
		if (fabs(nn-nstep/2) <= 2*k+1) {
			ipt++;
			for (int i=0;i<n;i++)
				fsave[ipt][i]=yend[i];
		}
		for (int i=0;i<n;i++) {
			double swap=ym[i]+h2*yend[i];
			ym[i]=yn[i];
			yn[i]=swap;
		}
		xnew += h;
		derivs(xnew,yn,yend);
	}
	if (nstep/2 <= 2*k+1) {
		ipt++;
		for (int i=0;i<n;i++)
			fsave[ipt][i]=yend[i];
	}
	for (int i=0;i<n;i++)
		yend[i]=0.5*(ym[i]+yn[i]+h*yend[i]);
}
template <class D>
void StepperBS<D>::polyextr(const int k, std::vector< std::vector<double> > &table, std::vector<double> &last) {
	int l=last.size();
	for (int j=k-1; j>0; j--)
		for (int i=0; i<l; i++)
			table[j-1][i]=table[j][i]+coeff[k][j]*(table[j][i]-table[j-1][i]);
	for (int i=0; i<l; i++)
		last[i]=table[0][i]+coeff[k][0]*(table[0][i]-last[i]);
}
template <class D>
void StepperBS<D>::step(const double htry,D &derivs) {
	const double STEPFAC1=0.65,STEPFAC2=0.94,STEPFAC3=0.02,STEPFAC4=4.0,
		KFAC1=0.8,KFAC2=0.9;
	static bool first_step=true,last_step=false;
	static bool forward,reject=false,prev_reject=false;
	int i,k;
	double fac,h,hnew,hopt_int,err;
	bool firstk;
	std::vector<double> hopt(IMAXX),work(IMAXX);
	std::vector<double> ysav(n),yseq(n);
	std::vector<double> ymid(n),scale(n);
	work[0]=0;
	h=htry;
	forward = h>0 ? true : false;
	for (i=0;i<n;i++) ysav[i]=y[i];
	if (h != hnext && !first_step) {
		last_step=true;
	}
	if (reject) {
		prev_reject=true;
		last_step=false;
	}
	reject=false;
	firstk=true;
	hnew=abs(h);
	interp_error:
	while (firstk || reject) {
		h = forward ? hnew : -hnew;
		firstk=false;
		reject=false;
		if (abs(h) <= abs(x)*EPS)
            std::cerr << "step size underflow in StepperBS" << std::endl;
		int ipt=-1;
		for (k=0; k<=k_targ+1;k++) {
			dy(ysav,h,k,yseq,ipt,derivs);
			if (k == 0)
				 y=yseq;
			else
				for (i=0;i<n;i++)
					table[k-1][i]=yseq[i];
			if (k != 0) {
				polyextr(k,table,y);
				err=0.0;
				for (i=0;i<n;i++) {
					scale[i]=atol+rtol*MAX(abs(ysav[i]),abs(y[i]));
					err+=SQR((y[i]-table[0][i])/scale[i]);
				}
				err=sqrt(err/n);
				double expo=1.0/(2*k+1);
				double facmin=pow(STEPFAC3,expo);
				if (err == 0.0)
					fac=1.0/facmin;
				else {
					fac=STEPFAC2/pow(err/STEPFAC1,expo);
					fac=MAX(facmin/STEPFAC4,MIN(1.0/facmin,fac));
				}
				hopt[k]=abs(h*fac);
				work[k]=cost[k]/hopt[k];
				if ((first_step || last_step) && err <= 1.0)
					break;
				if (k == k_targ-1 && !prev_reject && !first_step && !last_step) {
					if (err <= 1.0)
						break;
					else if (err>Q(nseq[k_targ]*nseq[k_targ+1]/(nseq[0]*nseq[0]))) {
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
					else if (err>SQR(nseq[k+1]/nseq[0])) {
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
		if (reject)
			prev_reject=true;
	}
	derivs(x+h,y,dydxnew);

	dydx=dydxnew;
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

#endif
