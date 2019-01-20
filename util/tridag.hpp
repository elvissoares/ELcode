#ifndef TRIDAG_HPP_
#define TRIDAG_HPP_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits> 

#include "nr.hpp"

void tridag(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &r, std::vector<double> &u)
//Resolve para um vetor u[0..n-1] o sistema linear tridiagonal dado pela equação (2.4.1 - NRC++). a[0..n-1], b[0..n-1], c[0..n-1] e r[0..n-1] são ventores input e não são modificados. b é a diagonal da matriz, enquanto que a[0] e c[n-1] são indefinidos e não são referenciados pela rotina que segue.
{
	int j,n=a.size();
	double bet;
	std::vector<double> gam(n);

	if (b[0] == 0.0) {
		std::cerr<< "ERROR: tridag: Error 1" << std::endl;
		exit(0);
	}
	u[0]=r[0]/(bet=b[0]);

	for (j=1;j<n;j++) {     //Decomposição e substituição para frente.
		gam[j]=c[j-1]/bet;
		bet=b[j]-a[j]*gam[j];
		if (bet == 0.0) {
			std::cerr<< "ERROR: tridag: Error 2" << std::endl;
			exit(0);
		}
		u[j]=(r[j]-a[j]*u[j-1])/bet;
	}
	for (j=(n-2);j>=0;j--)
		u[j] -= gam[j+1]*u[j+1];
}

void tridag_symm(std::vector<double> &c, std::vector<double> &b, std::vector<double> &r, std::vector<double> &u)
//Resolve para um vetor u[0..n-1] o sistema linear tridiagonal simétrico dado pela equação (2.4.1 - NRC++). c[0..n-1], b[0..n-1] e r[0..n-1] são ventores input e não são modificados. b é a diagonal da matriz, enquanto que c[n-1] é indefinido e não é referenciado pela rotina que segue.
{
	int n = b.size();
	std::vector<double> a(n);

	for (unsigned int j=1;j<n;j++)  //Copia o vetor c para a, de forma correta.
		a[j] = c[j-1];

	tridag(a,b,c,r,u);

	a.clear();
}

void cyclic(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c,const double alpha, const double beta, std::vector<double> &r, std::vector<double> &x)
{
	int i,n=a.size();
	double fact,gamma;
	if (n <= 2) throw("n too small in cyclic");
	std::vector<double> bb(n),u(n),z(n);
	gamma = -b[0];
	bb[0]=b[0]-gamma;
	bb[n-1]=b[n-1]-alpha*beta/gamma;
	for (i=1;i<n-1;i++) bb[i]=b[i];
	tridag(a,bb,c,r,x);
	u[0]=gamma;
	u[n-1]=alpha;
	for (i=1;i<n-1;i++) u[i]=0.0;
	tridag(a,bb,c,u,z);
	fact=(x[0]+beta*x[n-1]/gamma)/
		(1.0+z[0]+beta*z[n-1]/gamma);
	for (i=0;i<n;i++) x[i] -= fact*z[i];
}

#endif
