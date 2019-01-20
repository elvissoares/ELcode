#include "fermions.hpp"
#include "../cpp_constants_CGS.hpp"

#include "../util/roots.hpp"

Fermions::Fermions(){
    using PhysConstants::h;
    using PhysConstants::c;
    using PhysConstants::pi;
    using PhysConstants::M_e;
    using PhysConstants::k_B;
    using PhysConstants::N_A;

    B_mu = 8.*pi*C(M_e)*C(c)/(3.*N_A*C(h));
    A = pi*pow(M_e,4.)*pow(c,5.)/(3.*C(h));
}

//=====================================================================
// Variável auxiliar \beta que normaliza a temperatura em unidades
// da massa de repouso do férmion.
//=====================================================================
double Fermions::Beta(const double& T){
    
    using PhysConstants::c;
    using PhysConstants::k_B;
    using PhysConstants::M_e;
    
    return (k_B*T)/(M_e*Q(c));
}

void Fermions::EvaluateIntegrals(const double& eta, const double& beta)
{
	F12_n = F(0.5,eta,beta) - F(0.5,-eta-2/beta,beta) ;
    F32_n = F(1.5,eta,beta) - F(1.5,-eta-2/beta,beta) ;
    dF12deta_n = dFdeta(0.5,eta,beta) + dFdeta(0.5,-eta-2/beta,beta);
    dF32deta_n = dFdeta(1.5,eta,beta) + dFdeta(1.5,-eta-2/beta,beta);
    dF12dbeta_n = dFdbeta(0.5,eta,beta) - dFdbeta(0.5,-eta-2/beta,beta) - (2/Q(beta))*dFdeta(0.5,-eta-2/beta,beta) ;
    dF32dbeta_n = dFdbeta(1.5,eta,beta) - dFdbeta(1.5,-eta-2/beta,beta) - (2/Q(beta))*dFdeta(1.5,-eta-2/beta,beta);

    F32 = F(1.5,eta,beta) + F(1.5,-eta-2/beta,beta) ;
    F52 = F(2.5,eta,beta) + F(2.5,-eta-2/beta,beta) ;
    dF32deta = dFdeta(1.5,eta,beta) - dFdeta(1.5,-eta-2/beta,beta);
    dF52deta = dFdeta(2.5,eta,beta) - dFdeta(2.5,-eta-2/beta,beta);
    dF32dbeta = dFdbeta(1.5,eta,beta) + dFdbeta(1.5,-eta-2/beta,beta) + (2/Q(beta))*dFdeta(1.5,-eta-2/beta,beta);
    dF52dbeta = dFdbeta(2.5,eta,beta) + dFdbeta(2.5,-eta-2/beta,beta) + (2/Q(beta))*dFdeta(2.5,-eta-2/beta,beta);
}

void Fermions::Evaluate(const double& rho_in, const double& T_in, const double &Ye_in){

    using PhysConstants::N_A;
    using PhysConstants::k_B;

    rho = rho_in;
    T = T_in;
    Ye = Ye_in;
    
    n = rho*N_A*Ye; // number of free electrons

    eta = Eta_Search(n,T); // eta & beta parameters
    beta = Beta(T); 

    EvaluateIntegrals(eta,beta);

    p = Pressure();
    u = EnergyDensity()/rho; //specific internal energy
    s = ((u*rho+p)/T - eta*k_B*n)/rho; //specific entropy
    dpdrho = dPdRho();
    dpdT = dPdT();

    dudrho = dUdRho()/rho - u/rho; 
    dudT = dUdT()/rho;
    dsdrho = -(1/Q(rho))*dpdT;
    dsdT = dudT/T;
}



//=====================================================================
// Densidade de número de elétrons em função de \beta e \eta
// n_e(\beta,\eta), auxiliar para cálculo de \eta
//=====================================================================

double Fermions::Pressure()
{ 
    return (16.*sqrt(2)*A)*pow(beta,2.5)*(F32 + 0.5*beta*F52);
}

//=====================================================================
// Densidade de energia u para elétrons em função da densidade de
// massa \rho e da temperatura T
//=====================================================================
double Fermions::EnergyDensity()
{
	return 24.*sqrt(2)*A*pow(beta,1.5)*(beta*(F32 + beta*F52) 
    			+ 2*(F(0.5,-eta-2/beta,beta) + beta*F(1.5,-eta-2/beta,beta)) );

}

struct root_Eta : public Fermions{
    double np, Tp; // n-point and T-point

    double operator()(double etaa){

    	using PhysConstants::N_A;

        double betaa = Beta(Tp);
        double f12, f32;

		f12 = F(0.5,etaa,betaa) - F(0.5,-etaa-2/betaa,betaa) ;
    	f32 = F(1.5,etaa,betaa) - F(1.5,-etaa-2/betaa,betaa) ;

        double n_calc = 3*sqrt(2)*N_A*B_mu*pow(betaa,1.5)*(f12 + betaa*f32);
    
        return (log(np)-log(n_calc));
    }
};

double Fermions::Eta_Search(const double& np, const double& Tp)
{   
    root_Eta func;

    func.np = np;
    func.Tp = Tp;

    double betaa = Beta(Tp);

    double etaa = zbrent(func,-1/betaa,1e12,1e-16);

    return etaa;
}

double Fermions::dPdRho()
{ 
    using PhysConstants::N_A;
    using PhysConstants::M_e;
    using PhysConstants::c;

    return Ye*N_A*(2/3.)*M_e*Q(c)*beta*(dF32deta + 0.5*beta*dF52deta)/(dF12deta_n + beta*dF32deta_n);
}

double Fermions::dPdT()
{
    using PhysConstants::N_A;
    using PhysConstants::k_B;
    using PhysConstants::M_e;
    using PhysConstants::c;

    double dPdbeta = 16*sqrt(2)*A*pow(beta,1.5)*(2.5*F32 + (7/4.)*beta*F52 + beta*dF32dbeta + 0.5*Q(beta)*dF52dbeta);
    double dPdeta = 16*sqrt(2)*A*pow(beta,2.5)*(dF32deta + 0.5*beta*dF52deta);

    double dndbeta = 3*sqrt(2)*N_A*(B_mu)*sqrt(beta)*( 1.5*F12_n + 2.5*beta*F32_n + beta*dF12dbeta_n + Q(beta)* dF32dbeta_n);
    double dndeta = 3*sqrt(2)*N_A*(B_mu)*pow(beta,1.5)*(dF12deta_n + beta*dF32deta_n);

    return (k_B/(M_e*Q(c)))*(dPdbeta- dPdeta*dndbeta/dndeta);
}

double Fermions::dUdRho()
{ 
    using PhysConstants::N_A;
    using PhysConstants::M_e;
    using PhysConstants::c;

    double dudeta = 24*sqrt(2)*A*pow(beta,1.5)*( beta*(dF32deta + beta*dF52deta)
           - 2*(dFdeta(0.5,-eta-2/beta,beta) + beta*dFdeta(1.5,-eta-2/beta,beta)) );
    double dndeta = 3*sqrt(2)*N_A*(B_mu)*pow(beta,1.5)*(dF12deta_n + beta*dF32deta_n);

    return Ye*N_A*dudeta/dndeta;
}

double Fermions::dUdT()
{
    using PhysConstants::N_A;
    using PhysConstants::k_B;
    using PhysConstants::M_e;
    using PhysConstants::c;
    
    double dUdbeta = 24*sqrt(2)*A*sqrt(beta)*( beta*(2.5*F32 + (3.5)*beta*F52 + beta*dF32dbeta + Q(beta)*dF52dbeta) + 2*(1.5*F(0.5,-eta-2/beta,beta) + 2.5*beta*F(1.5,-eta-2/beta,beta) + 2*dFdeta(0.5,-eta-2/beta,beta)/beta + beta*dFdbeta(0.5,-eta-2/beta,beta) + 2*beta*dFdeta(1.5,-eta-2/beta,beta) + Q(beta)*dFdbeta(1.5,-eta-2/beta,beta)) );
    
    double dUdeta = 24*sqrt(2)*A*pow(beta,1.5)*( beta*(dF32deta + beta*dF52deta)
                                                - 2*(dFdeta(0.5,-eta-2/beta,beta) + beta*dFdeta(1.5,-eta-2/beta,beta)) );
    
    double dndbeta = 3*sqrt(2)*N_A*(B_mu)*sqrt(beta)*( 1.5*F12_n + 2.5*beta*F32_n + beta*dF12dbeta_n + Q(beta)* dF32dbeta_n);
    double dndeta = 3*sqrt(2)*N_A*(B_mu)*pow(beta,1.5)*(dF12deta_n + beta*dF32deta_n);
    
    return (k_B/(M_e*Q(c)))*(dUdbeta- dUdeta*dndbeta/dndeta);
}

double Fermions::Free(const double &rhoin, const double &Tin)
{
    using PhysConstants::N_A;
	using PhysConstants::k_B;

    double nn = Ye*rhoin*N_A; // number of free electrons

    double etaa = Eta_Search(nn,Tin); // eta & beta parameters
    double betaa = Beta(Tin); 

    double f32 = F(1.5,etaa,betaa) + F(1.5,-etaa-2/betaa,betaa) ;
    double f52 = F(2.5,etaa,betaa) + F(2.5,-etaa-2/betaa,betaa) ;

    double pp = (16.*sqrt(2)*A)*pow(betaa,2.5)*(f32 + 0.5*betaa*f52);
    double uu = 24.*sqrt(2)*A*pow(betaa,1.5)*(betaa*(f32 + betaa*f52) 
    			+ 2*(F(0.5,-etaa-2/betaa,betaa) + betaa*F(1.5,-etaa-2/betaa,betaa)) );

    return ((uu+pp)/Tin - etaa*k_B*nn)/rhoin;
}

void Fermions::Evaluate_FreeEnergyDerivatives()
{
	const int ntab=40;
	const double con=1.4, con2=(con*con);
	const double big=std::numeric_limits<double>::max();
	const double safe=2.0;

	double errta, errtb, errtc ,fac, h= 0.05*rho, k = 0.05*T,
			ansa, ansb, ansc;

	std::vector< std::vector<double> > a, b ,c;
	a.resize(ntab);
	b.resize(ntab);
	c.resize(ntab);
	for (unsigned int i=0;i<ntab;i++) {
		a[i].resize(ntab);
		b[i].resize(ntab);
		c[i].resize(ntab);
	}

	double pt0 = Free(rho,T);
	double pt1 = Free(rho,T+k);
	double pt2 = Free(rho,T-k);
	double pt3 = Free(rho-h,T);
	double pt4 = Free(rho+h,T);
	double pt5 = Free(rho-h,T-k);
	double pt6 = Free(rho-h,T+k);
	double pt7 = Free(rho+h,T-k);
	double pt8 = Free(rho+h,T+k);

    a[0][0]=((pt7-2*pt4+pt8)-(pt5-2*pt3+pt6))/(2*h*k*k); //d3FdT2dRho
    b[0][0]=((pt6-2*pt1+pt8)-(pt5-2*pt2+pt7))/(2*h*h*k); //d3FdRho2dT
    c[0][0]= ( (pt5 - 2*pt3 + pt6) - 2*(pt2 - 2*pt0 + pt1) + (pt7 - 2*pt4 + pt8))/(h*h*k*k);
	double err=big;

	for (unsigned int i=1;i<ntab;i++) {
		h /= con;
		k /= con;

		pt0 = Free(rho,T);
		pt1 = Free(rho,T+k);
		pt2 = Free(rho,T-k);
		pt3 = Free(rho-h,T);
		pt4 = Free(rho+h,T);
		pt5 = Free(rho-h,T-k);
		pt6 = Free(rho-h,T+k);
		pt7 = Free(rho+h,T-k);
		pt8 = Free(rho+h,T+k);

		a[0][i]=((pt7-2*pt4+pt8)-(pt5-2*pt3+pt6))/(2*h*k*k); //d3FdT2dRho
		b[0][i]=((pt6-2*pt1+pt8)-(pt5-2*pt2+pt7))/(2*h*h*k); //d3FdRho2dT
		c[0][i]= ( (pt5 - 2*pt3 + pt6) - 2*(pt2 - 2*pt0 + pt1) + (pt7 - 2*pt4 + pt8))/(h*h*k*k); 
		fac=con2;

		for (unsigned int j=1;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			b[j][i]=(b[j-1][i]*fac-b[j-1][i-1])/(fac-1.0);
			c[j][i]=(c[j-1][i]*fac-c[j-1][i-1])/(fac-1.0);
			fac=con2*fac;

			errta=MAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
			errtb=MAX(fabs(b[j][i]-b[j-1][i]),fabs(b[j][i]-b[j-1][i-1]));
			errtc=MAX(fabs(c[j][i]-c[j-1][i]),fabs(c[j][i]-c[j-1][i-1]));

			if (errta <= err && errtb <= err && errtc <= err) {
				err=MAX(errta,MAX(errtb,errtc));
				ansa=a[j][i];
				ansb=b[j][i];
				ansc=c[j][i];
			}
		}

		if (fabs(a[i][i]-a[i-1][i-1]) >= safe*err && fabs(b[i][i]-b[i-1][i-1]) >= safe*err &&
			fabs(c[i][i]-c[i-1][i-1]) >= safe*err ) break;
	}

	d3fdT2drho = ansa;
	d3fdrho2dT= ansb;
	d4fdrho2dT2 = ansc;
}

void InterpolatedFermions::Evaluate(const double &Rhoin, const double &Tin, const double &Ye)
{
	double dT = T[1] - T[0];
	double drho = Rho[1]-Rho[0];
	//Find the table locations (the table have to be log spaced)
	int jat = int((log10(Tin)-T[0])/dT);
	double din = Ye*Rhoin;
	int iat = int((log10(din)-Rho[0])/drho);

	if (log10(Tin) > T.back() || log10(din) > Rho.back()) {
		std::cerr << "Fermions: out range of table" << std::endl;

		Fermions fermions;
		fermions.Evaluate(Rhoin,Tin,Ye);

		p = fermions.p;
		dpdT = fermions.dpdT;
		dpdrho = fermions.dpdrho;

		s = fermions.s;
		dsdT = fermions.dsdT;
    	dsdrho = fermions.dsdrho;

		u = fermions.u;
		dudT = fermions.dudT;
		dudrho = fermions.dudrho;

		return;
	}

	//various differences
	double dt = pow(10,T[jat+1]) - pow(10,T[jat]);
	double dt2 = dt*dt;
	double dd = pow(10,Rho[iat+1]) - pow(10,Rho[iat]);
	double dd2 = dd * dd;
	double xt = MAX((Tin-pow(10,T[jat]))/dt,0.0);
	double xd = MAX((din-pow(10,Rho[iat]))/dd,0.0);
	double mxt = 1.0 - xt;
	double mxd = 1.0 - xd;

	//evaluate the basis function
	double si0t = psi0(xt);
	double si1t = psi1(xt)*dt;
	double si2t = psi2(xt)*dt2;

	double si0mt = psi0(mxt);
	double si1mt = -psi1(mxt)*dt;
	double si2mt = psi2(mxt)*dt2;

	double si0d = psi0(xd);
	double si1d = psi1(xd)*dd;
	double si2d = psi2(xd)*dd2;

	double si0md = psi0(mxd);
	double si1md = -psi1(mxd)*dd;
	double si2md = psi2(mxd)*dd2;

	//and their first derivatives
	double dsi0t = dpsi0(xt)/dt;
	double dsi1t = dpsi1(xt);
	double dsi2t = dpsi2(xt)*dt;

	double dsi0mt = -dpsi0(mxt)/dt;
	double dsi1mt = dpsi1(mxt);
	double dsi2mt = -dpsi2(mxt)*dt;

	double dsi0d = dpsi0(xd)/dd;
	double dsi1d = dpsi1(xd);
	double dsi2d = dpsi2(xd)*dd;

	double dsi0md = -dpsi0(mxd)/dd;
	double dsi1md = dpsi1(mxd);
	double dsi2md = -dpsi2(mxd)*dd;

	//and their second derivatives
	double ddsi0t = ddpsi0(xt)/dt2;
	double ddsi1t = ddpsi1(xt)/dt;
	double ddsi2t = ddpsi2(xt);

	double ddsi0mt = ddpsi0(mxt)/dt2;
	double ddsi1mt = -ddpsi1(mxt)/dt;
	double ddsi2mt = ddpsi2(mxt);

	double ddsi0d = ddpsi0(xd)/dd2;
	double ddsi1d = ddpsi1(xd)/dd;
	double ddsi2d = ddpsi2(xd);

	double ddsi0md = ddpsi0(mxd)/dd2;
	double ddsi1md = -ddpsi1(mxd)/dd;
	double ddsi2md = ddpsi2(mxd);

	//the free energy
	double Free = herm5(iat,jat,si0t,si1t,si2t,si0mt,si1mt,si2mt,
					si0d,si1d,si2d,si0md,si1md,si2md);

	//the derivative of the free energy with density
	double dFdrho = herm5(iat,jat,si0t,si1t,si2t,si0mt,si1mt,si2mt,
						dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md);

	//the derivative of the free energy with temperature
	double dFdT = herm5(iat,jat,dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt,
					si0d,si1d,si2d,si0md,si1md,si2md);

	//the second derivative of the free energy with density^{-2}
	double d2Fdrho2 = herm5(iat,jat,si0t,si1t,si2t,si0mt,si1mt,si2mt,
						ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md);

	//the second derivative of the free energy with temperature^{-2}
	double d2FdT2 = herm5(iat,jat,ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt,
					si0d,si1d,si2d,si0md,si1md,si2md);

	//the second derivative of the free energy with density and temperasture
	double d2FdrhodT = herm5(iat,jat,dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt,
						dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md);

	//output: pressure, specific energy, specific entropy and their partial derivatives
	p = Q(din)*dFdrho;
	dpdT = Q(din)*d2FdrhodT;
	dpdrho = Ye*(Q(din)*d2Fdrho2 + 2.*din*dFdrho);

	s = - dFdT * Ye;
	dsdT = -d2FdT2 * Ye;
    dsdrho = -d2FdrhodT * Q(Ye);

	u = Ye * Free + Tin *s;
	dudT = Tin * dsdT;
	dudrho = Q(Ye) * dFdrho + Tin *dsdrho;
}

// Quintic Hermite basis functin Psi0
double InterpolatedFermions::psi0(const double &z)
{
	return (-6.*pow(z,5) + 15.*pow(z,4) - 10*pow(z,3) + 1.);
}

double InterpolatedFermions::dpsi0(const double &z)
{
	return (-30.*pow(z,4) + 60.*pow(z,3) - 30.*pow(z,2));
}

double InterpolatedFermions::ddpsi0(const double &z)
{
	return (-120.*pow(z,3) + 180.*pow(z,2) - 60.*z);
}

// Quintic Hermite basis functin Psi1
double InterpolatedFermions::psi1(const double &z)
{
	return (-3.*pow(z,5) + 8*pow(z,4) - 6.*pow(z,3) + z);
}

double InterpolatedFermions::dpsi1(const double &z)
{
	return (-15.*pow(z,4) + 32.*pow(z,3) - 18.*pow(z,2) + 1.);
}

double InterpolatedFermions::ddpsi1(const double &z)
{
	return (-60.*pow(z,3) + 96.*pow(z,2) - 36.*z);
}

// Quintic Hermite basis functin Psi2
double InterpolatedFermions::psi2(const double &z)
{
	return 0.5*(-pow(z,5) +3.*pow(z,4) -3.*pow(z,3) +  pow(z,2));
}

double InterpolatedFermions::dpsi2(const double &z)
{
	return 0.5*(-5.*pow(z,4) + 12.*pow(z,3) - 9.*pow(z,2) + 2.*z);
}

double InterpolatedFermions::ddpsi2(const double &z)
{
	return 0.5*(-20.*pow(z,3) + 36.*pow(z,2) - 18.*z + 2.);
}

double InterpolatedFermions::herm5(const int &i, const int &j, const double &w0t,
	const double &w1t, const double &w2t, const double &w0mt, const double &w1mt, 
	const double &w2mt, const double &w0d, const double &w1d, const double &w2d,
	const double &w0md, const double &w1md, const double &w2md)
{
	return (f[i][j]*w0d*w0t + f[i+1][j]*w0md*w0t + f[i][j+1]*w0d*w0mt + f[i+1][j+1]*w0md*w0mt
	+ ft[i][j]*w0d*w1t + ft[i+1][j]*w0md*w1t + ft[i][j+1]*w0d*w1mt + ft[i+1][j+1]*w0md*w1mt
	+ fd[i][j]*w1d*w0t + fd[i+1][j]*w1md*w0t + fd[i][j+1]*w1d*w0mt + fd[i+1][j+1]*w1md*w0mt
	+ ftt[i][j]*w0d*w2t + ftt[i+1][j]*w0md*w2t + ftt[i][j+1]*w0d*w2mt + ftt[i+1][j+1]*w0md*w2mt
	+ fdd[i][j]*w2d*w0t + fdd[i+1][j]*w2md*w0t + fdd[i][j+1]*w2d*w0mt + fdd[i+1][j+1]*w2md*w0mt
    + fdt[i][j]*w1d*w1t + fdt[i+1][j]*w1md*w1t + fdt[i][j+1]*w1d*w1mt + fdt[i+1][j+1]*w1md*w1mt
	+ fdtt[i][j]*w1d*w2t + fdtt[i+1][j]*w1md*w2t + fdtt[i][j+1]*w1d*w2mt + fdtt[i+1][j+1]*w1md*w2mt
	+ fddt[i][j]*w2d*w1t + fddt[i+1][j]*w2md*w1t + fddt[i][j+1]*w2d*w1mt + fddt[i+1][j+1]*w2md*w1mt
	+ fddtt[i][j]*w2d*w2t + fddtt[i+1][j]*w2md*w2t + fddtt[i][j+1]*w2d*w2mt + fddtt[i+1][j+1]*w2md*w2mt);
}

void InterpolatedFermions::Read_Table()
{

	std::ifstream InFile("input/free_energy-pairs.dat");

	unsigned int N, M;

	std::string Line;
    
    if (InFile.is_open()){
		
		getline(InFile, Line);
		
		getline(InFile, Line);
		std::string post = Line.substr(Line.find("=")+1);
		N = atof(post.c_str());
		
		getline(InFile, Line);
		post = Line.substr(Line.find("=")+1);
		M = atof(post.c_str());

		getline(InFile, Line);
    }
    
    else {
       std::cerr << "InterpolatedFermions: The input file doesnt exist!";
        exit(1);
    }

	Rho.resize(N); T.resize(M);
    
    // Set up sizes. (HEIGHT x WIDTH)
    f.resize(N); ft.resize(N); fd.resize(N);
    ftt.resize(N); fdd.resize(N); fdt.resize(N); 
    fddt.resize(N); fdtt.resize(N); fddtt.resize(N);
        
	for (unsigned int i = 0; i < N; ++i){
		f[i].resize(M); ft[i].resize(M); fd[i].resize(M);
    	ftt[i].resize(M); fdd[i].resize(M); fdt[i].resize(M); 
    	fddt[i].resize(M); fdtt[i].resize(M); fddtt[i].resize(M);
	}
	
	unsigned int i = 0, j = 0;
    
    while(j < M){
		InFile >> Rho[i] >> T[j] >> f[i][j] >> fd[i][j] >> ft[i][j] 
		>> fdd[i][j] >> ftt[i][j] >> fdt[i][j] >> fddt[i][j] >> fdtt[i][j] >> fddtt[i][j];

		// Rhomin = MIN(Rhomin,Rho[i]);
		// Rhomax = MAX(Rhomax,Rho[i]);
		// Tmin = MIN(Tmin,T[i]);
		// Tmax = MAX(Tmax,T[i]);

		++i;
		if (i==N) {
            i=0;	++j;
            getline(InFile, Line);
        }
        
    }

}
