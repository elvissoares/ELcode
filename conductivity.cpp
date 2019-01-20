#include "conductivity.hpp"
#include "cpp_constants_CGS.hpp"
#include "eos.hpp"
#include "util/nr.hpp"

double E1(const double &x){
	//Computa a integral de Fermi-Dirac generalizada F_k(eta,theta) onde k > -1 e \theta >= 0. A acurácia é aproximada-
	// mente o quadrado do parâmetro EPS. NMAX limita o número total de passos da quadratura.
	const double EPS = 1.e-8;
	const double TINY = std::numeric_limits<double>::epsilon();
	const int NMAX = 16;
	double a, aa, b, bb, hmax, olds, sum;

	a=0.0;		//Atribua limites para o mapeamento x=exp(t-e^{-t})
	b=exp(-x);
	ExponentialIntegral expint;
	Trapzd<ExponentialIntegral> s(expint,a,b);
	for (unsigned int i=1;i<=NMAX;i++){
		sum=s.next();
		if (i>3){ 		//Teste para convergência
			if(fabs(sum-olds) <= EPS*fabs(olds) + TINY)
				return sum;
		}
		olds=sum; //Valor salvo para próximo teste de convergência
	}
}


int hunt(const double x, const std::vector<double> xx){

	unsigned int n = xx.size();

	bool ascnd = (xx[n-1] >= xx[0]); //true if ascending order, false otherwise

	int jlo=n/2+1, jhi, jm, inc=1;

	if (jlo < 0 || jlo > n-1){
		jlo = 0.;
		jhi = n-1;
	}

	else {
		if (x>= xx[jlo] == ascnd){
			for (;;){
				jhi = jlo + inc;
				if (jhi >= n-1) { jhi = n-1; break;}
				else if (x < xx[jhi] == ascnd) break;
				else {
					jlo = jhi;
					inc += inc;
				}
			}
		}
		else{
			jhi = jlo;
			for (;;){
				jlo = jlo - inc;
				if (jlo <= 0) { jlo = 0; break;}
				else if (x >= xx[jlo] == ascnd) break;
				else {
					jhi = jlo;
					inc += inc;
				}
			}
		}
	}
	while (jhi-jlo > 1){
		jm = (jhi+jlo) >> 1;
		if (x >= xx[jm] == ascnd)
			jlo=jm;
		else jhi = jm;
	}

	return jlo;

}

void bcucof(const std::vector<double> &y, const std::vector<double> &y1, const std::vector<double> &y2, const std::vector<double> &y12, const double &d1, const double &d2, std::vector< std::vector<double> > &c)
{

	static int wt_d[16*16] = 
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

	double xx, d1d2 = d1*d2;
	std::vector<double> cl(16), x(16);
	std::vector< std::vector<double> > wt(16,std::vector<double>(16));

	for (unsigned int i = 0; i<16; i++){ 
		for (unsigned int j = 0; j<16; j++)
			wt[i][j] = wt_d[i*16+j];	
	}

	for (unsigned int i = 0; i<4; i++){ //Carregue um vetor temporário
		x[i] = y[i];
		x[i+4] = y1[i]*d1;
		x[i+8]= y2[i]*d2;
		x[i+12] = y12[i]*d1d2;
	}

	for (unsigned int i = 0; i<16; i++){
		xx = 0.0;
		for (unsigned int k = 0; k<16; k++) xx += wt[i][k] * x[k];
		cl[i] = xx;
	}

	int l = 0;
	for (unsigned int i = 0; i<4; i++)
		for (unsigned int j = 0; j<4; j++) c[i][j] = cl[l++];
}

void bcuint(const std::vector<double> &y, const std::vector<double> &y1, const std::vector<double> &y2, const std::vector<double> &y12, const double &x1lo, const double &x1hi, const double &x2lo, const double &x2hi, const double &x1, const double &x2, double &ansy, double &ansy1, double &ansy2)
{
	double t, u, d1 = x1hi-x1lo, d2=x2hi-x2lo;
	std::vector< std::vector<double> > c(4,std::vector<double>(4));

	if (x1lo == x1hi || x2lo == x2hi){
		std::cerr << "Bad input in routine bcuint" << std::endl;
		exit(1);
	}

	bcucof(y,y1,y2,y12,d1,d2,c);

	t = (x1-x1lo)/d1;
	u = (x2-x2lo)/d2;

	ansy = ansy2 = ansy1 = 0.0;

	for (unsigned int i = 0; i<4; i++){
		for (unsigned int j = 0; j<4; j++)
			ansy+=c[i][j]*pow(t,i)*pow(u,j);
	}

	ansy1 /= d1;
	ansy2 /= d2;
}


void ThermalConductivity::Read_Table()
{

	std::ifstream InFile("input/condall06.dat");

	unsigned int MAXT=19, MAXR=64, MAXZ=15;

	RHO.resize(MAXR); T.resize(MAXT); Z.resize(MAXZ);
    
    // Set up sizes. (HEIGHT x WIDTH)
    KAP.resize(MAXZ);
        
	for (unsigned int i = 0; i < MAXZ; ++i){
		KAP[i].resize(MAXR);
		for (unsigned int j = 0; j < MAXR; ++j)
			KAP[i][j].resize(MAXT);
	}

	std::string Line;
    
    if (InFile.is_open()){
		
		getline(InFile, Line);// Skip the first line
    }
    
    else {
       std::cerr << "Conductivity: The input file doesnt exist!";
        exit(1);
    }
    
    for (unsigned int i = 0; i < MAXZ; ++i){
		InFile >> Z[i];
		for (unsigned int j = 0; j < MAXT; ++j){
		 	InFile >> T[j];
		}

		for (unsigned int k = 0; k < MAXR; ++k){
			InFile >> RHO[k];
			for (unsigned int j = 0; j < MAXT; ++j){
		 		InFile >> KAP[i][k][j];
			}
		}
	}

}

void ThermalConductivity::Num_Derivative(std::vector<double> &y, std::vector<double> &y1, std::vector<double> &y2, std::vector<double> &y12, const int &iZ, const int &iR, const int &iT)
{
	y[0] = KAP[iZ][iR][iT];
	y[1] = KAP[iZ][iR+1][iT];
	y[2] = KAP[iZ][iR+1][iT+1];
	y[3] = KAP[iZ][iR][iT+1];

	y1[0] = (KAP[iZ][iR+1][iT]- KAP[iZ][iR-1][iT])/ (RHO[iR+1]-RHO[iR-1]);
	y1[1] = (KAP[iZ][iR+2][iT]- KAP[iZ][iR][iT])/ (RHO[iR+2]-RHO[iR]);
	y1[2] = (KAP[iZ][iR+2][iT+1]- KAP[iZ][iR][iT+1])/ (RHO[iR+2]-RHO[iR]);
	y1[3] = (KAP[iZ][iR+1][iT+1]- KAP[iZ][iR-1][iT+1])/ (RHO[iR+1]-RHO[iR-1]);

	y2[0] = (KAP[iZ][iR][iT+1]- KAP[iZ][iR][iT-1])/ (T[iT+1]-T[iT-1]);
	y2[1] = (KAP[iZ][iR+1][iT+1]- KAP[iZ][iR+1][iT-1])/ (T[iT+1]-T[iT-1]);
	y2[2] = (KAP[iZ][iR+1][iT+2]- KAP[iZ][iR+1][iT])/ (T[iT+2]-T[iT]);
	y2[3] = (KAP[iZ][iR][iT+2]- KAP[iZ][iR][iT])/ (T[iT+2]-T[iT]);

	y12[0] = (KAP[iZ][iR+1][iT+1] - KAP[iZ][iR+1][iT-1] - KAP[iZ][iR-1][iT+1] + KAP[iZ][iR-1][iT-1])/((RHO[iR+1]-RHO[iR-1])*(T[iT+1]-T[iT-1]));
	y12[1] = (KAP[iZ][iR+2][iT+1] - KAP[iZ][iR+2][iT-1] - KAP[iZ][iR][iT+1] + KAP[iZ][iR][iT-1])/((RHO[iR+2]-RHO[iR])*(T[iT+1]-T[iT-1]));
	y12[2] = (KAP[iZ][iR+2][iT+2] - KAP[iZ][iR+2][iT] - KAP[iZ][iR][iT+2] + KAP[iZ][iR][iT])/((RHO[iR+2]-RHO[iR])*(T[iT+2]-T[iT]));
	y12[3] = (KAP[iZ][iR+1][iT+2] - KAP[iZ][iR+1][iT] - KAP[iZ][iR-1][iT+2] + KAP[iZ][iR-1][iT])/((RHO[iR+1]-RHO[iR-1])*(T[iT+2]-T[iT]));

}


double ThermalConductivity::Evaluate_ThermalConductivity(const double &Rhoin, const double &Tin, const std::vector<double> &Xin)
{
	//Rhoin = Log10(Rho) and Tin = Log10(T)
    static double Ze[7] = {2,6,8,10,12,14,28};
    static double A[7] = {4,12,16,20,24,28,56};
    double Zin = Ze[0]*Xin[0] + Ze[1]*Xin[1] + Ze[2]*Xin[2]+ Ze[3]*Xin[3] + Ze[4]*Xin[4]+ Ze[5]*Xin[5] + Ze[6]*Xin[6];
    double Ye = Ze[0]*Xin[0]/A[0] + Ze[1]*Xin[1]/A[1] + Ze[2]*Xin[2]/A[2]  + Ze[3]*Xin[3]/A[3]  + Ze[4]*Xin[4]/A[4] + Ze[5]*Xin[5]/A[5] + Ze[6]*Xin[6]/A[6] ; 
    double rhomu_e = Rhoin + log10(Ye) ;
    
	if (Tin > T.back() || Rhoin > RHO.back() ) {

		double Tdeg = (0.36*rhomu_e +7.34) ;

		if ( Tin > 1.1*Tdeg) return Nondegenerateconductivity(Rhoin,Tin,Xin);
		else if (Tin < 0.9*Tdeg ) return Degenerateconductivity(Rhoin,Tin,Xin);

		else {
			double a = Nondegenerateconductivity(Rhoin,Tin,Xin);
			double b = Degenerateconductivity(Rhoin,Tin,Xin);
			return 5.0*(a*(Tin - 0.9*Tdeg ) + b*(1.1*Tdeg -Tin))/Tdeg;
		}

	}

	int iZ = hunt(Zin,Z);
    int iR = hunt(Rhoin,RHO);
    int iT = hunt(Tin,T);

    std::vector<double> y(4), y1(4), y2(4), y12(4);

    // Z low = Z[iZ]
    Num_Derivative(y,y1,y2,y12,iZ,iR,iT);

    double kappaZlo, kappaZlo1, kappaZlo2;
    bcuint(y,y1,y2,y12,RHO[iR],RHO[iR+1],T[iT],T[iT+1],Rhoin,Tin,kappaZlo,kappaZlo1,kappaZlo2);

    // Z high = Z[iZ+1]
    Num_Derivative(y,y1,y2,y12,iZ+1,iR,iT);

    double kappaZhi, kappaZhi1, kappaZhi2;
    bcuint(y,y1,y2,y12,RHO[iR],RHO[iR+1],T[iT],T[iT+1],Rhoin,Tin,kappaZhi,kappaZhi1,kappaZhi2);

    //Por último realizamos uma interpolação linear em Z
    return  (kappaZlo + ((Zin - Z[iZ])/(Z[iZ+1]-Z[iZ]))*(kappaZhi - kappaZlo));
}

double ThermalConductivity::Nondegenerateconductivity(const double &logRho, const double &logT, const std::vector<double> &X)
{
	// This conductivity corresponds to the Appendix A in Iben, I., J. (1975). Thermal pulses; p-capture, alpha-capture, s-process nucleosynthesis; and convective mixing in a star of intermediate mass. The Astrophysical Journal, 196(9), 525. http://doi.org/10.1086/153433
	using PhysConstants::a;
	using PhysConstants::c;

	double Ze[7] = {2,6,8,10,12,14,28};
	static double A[7] = {4,12,16,20,24,28,56};

	double Ye = Ze[0]*X[0]/A[0] + Ze[1]*X[1]/A[1] + Ze[2]*X[2]/A[2]  + Ze[3]*X[3]/A[3]  + Ze[4]*X[4]/A[4] + Ze[5]*X[5]/A[5] + Ze[6]*X[6]/A[6] ; 
	double mu_e = 1.0/Ye;

	double T = pow(10.0,logT);
	double Rho = pow(10.0,logRho);
	double T6 = T/1.0e6;

	double delta = Rho*pow(T6,-1.5)/mu_e;
	double logdelta = log10(delta);

	//Calculation of the logtheta parameter
	double eta0 = pow(10.0,-0.52255+2*logdelta/3.0);
	double logtheta;

	if (logdelta <= 0.645) logtheta = -3.2862 + log10(delta*(1+0.024417*delta));
	else if (logdelta > 0.645 && logdelta < 2.0) logtheta = -3.29243 + log10(delta*(1+0.02804*delta));

	else if (logdelta > 2.5){
		logtheta = -4.80946 + log10(Q(delta)*(1+9.376/Q(eta0)));
	}
	else {
		double a1 = -3.29243 + log10(delta*(1+0.02804*delta));
		double b1 = -4.80946 + log10(Q(delta)*(1+9.376/Q(eta0)));
		logtheta = 2*a1*(2.5-logdelta) + 2*b1*(logdelta-2.0);
	}

	//Calculation of the (P/nkT) parameter
	double PovernkT;
	if (logdelta < 1.5) PovernkT = 1+0.021876*delta;
	else if (logdelta > 2.0) PovernkT = 0.4*eta0*(1+4.1124/Q(eta0));
	else {
		double a2 = log10(1+0.021876*delta);
		double b2 = log10(0.4*eta0*(1+4.1124/Q(eta0)));
		PovernkT = pow(10.0, 2*a2*(2.0-logdelta)+2*b2*(logdelta-1.5));
	}

	//Calculation of the n^{-1}\frac{\partial n_e}{\partial \eta} parameter
	double ndndeta;
	if (delta < 40.0) ndndeta = 1.0 - 0.01*delta*(2.8966-0.034838*delta);
	else ndndeta = (1.5/eta0)*(1.0-0.8225/Q(eta0));

	double XZ2A = Q(Ze[0])*X[0]/A[0] + Q(Ze[1])*X[1]/A[1]+ Q(Ze[2])*X[2]/A[2] + Q(Ze[3])*X[3]/A[3] + Q(Ze[4])*X[4]/A[4] + Q(Ze[5])*X[5]/A[5]+ Q(Ze[6])*X[6]/A[6]; 

	double lambda2Rsqr = (9.24735e-3)*(delta/PovernkT)*pow(T6,-0.5)*(mu_e*XZ2A+ndndeta);

	double alpha = 0.5*log10(lambda2Rsqr);

	double logthetaY;
	if (alpha <= -3.0) logthetaY = 0.937-0.111*alpha;
	else if (alpha > 0.0) logthetaY = 0.24-0.6*alpha;
	else logthetaY = 0.24-alpha*(0.55+0.0689*alpha);

	double logthetaC;
	if (alpha <= -2.5) logthetaC = 1.27-0.1*alpha;
	else if (alpha > 0.5) logthetaC = 0.843-0.785*alpha;
	else logthetaC = 0.727-alpha*(0.511+0.0778*alpha);

	double Zc = (A[1]*X[1]+A[2]*X[2]+A[3]*X[3]+A[4]*X[4]+A[5]*X[5]+A[6]*X[6])/12.0;

	double f = pow(10.0,logtheta);
	double thetaY = pow(10.0,logthetaY);
	double thetaC = pow(10.0,logthetaC);

	double opac = (X[0]*thetaY+Zc*thetaC)/(T6*f);

	return log10(4*a*c*C(T)/(3*Rho*opac));
}

double ThermalConductivity::Degenerateconductivity(const double &logRho, const double &logT, const std::vector<double> &X)
{

	//Contribution due to ee collisions, the reference is Shternin, P. S., & Yakovlev, D. G. (2006). Electron thermal conductivity owing to collisions between degenerate electrons. Physical Review D, 74(4), 43004. http://doi.org/10.1103/PhysRevD.74.043004
	using PhysConstants::hbar;
	using PhysConstants::pi;
	using PhysConstants::N_A;
	using PhysConstants::c;
	using PhysConstants::e_c;
	using PhysConstants::k_B;
    using PhysConstants::M_e;

	double Rho = pow(10.0,logRho);
	double T = pow(10.0,logT);
	double T6 = T/1.0e6;
	double Ze[7] = {2,6,8,10,12,14,28};
	double Ai[7] = {4,12,16,20,24,28,56};
	double Ye = Ze[0]*X[0]/Ai[0] + Ze[1]*X[1]/Ai[1] + Ze[2]*X[2]/Ai[2]  + Ze[3]*X[3]/Ai[3]  + Ze[4]*X[4]/Ai[4] + Ze[5]*X[5]/Ai[5] + Ze[6]*X[6]/Ai[6] ; 
	double Yi = X[0]/Ai[0] + X[1]/Ai[1] + X[2]/Ai[2] + X[3]/Ai[3] + X[4]/Ai[4] + X[5]/Ai[5] + X[6]/Ai[6] ; 

	double n_e = Rho*N_A*Ye;
	double n_i = Rho*N_A*Yi;
	double p_F = hbar*pow(3*Q(pi)*n_e,1/3.0);
	double k_F = p_F/hbar;
    
    double x_r = p_F/(M_e*c);
    double gamma_r = sqrt(1+Q(x_r));

    double m_estar = M_e*gamma_r;
    double v_F = p_F/m_estar;
    double T_F = (m_estar-M_e)*Q(c)/k_B;
    
	double omega_pe = sqrt(4*pi*Q(e_c)*n_e/m_estar);

	//Auxiliary variables from Yakovlev
	double beta =  x_r/gamma_r;
	double T_p = hbar*omega_pe/k_B;
	double y = sqrt(3.0)*T_p/T;

	double A = 12.2 + 25.2*C(beta);
    double C = 8.0/150+0.05714*pow(beta,4.0);
    double B = A*exp((0.123636+0.016234*Q(beta))/C);
	double D = 0.1558*pow(y,1-0.75*beta);

	double I;
//	if ( 1/y < 0.005) I = (10.0/63-(8.0/315)/(1+0.0435*y))*(128.56/(37.1*y+10.83*Q(y)+C(y)))/beta + C(beta)*(2.404/B + (C-2.404/B)/(1+0.1*y*beta))*(B/(A*y*beta+Q(y*beta))) + (beta/(1.0+D))*(C+18.52*Q(beta)*D/B)*(B/(A*y+10.83*Q(y*beta)+pow(y*beta,8/3.0)));
    
    I = (10.0/63-(8.0/315)/(1+0.0435*y))*log(1+128.56/(37.1*y+10.83*Q(y)+C(y)))/beta + C(beta)*(2.404/B + (C-2.404/B)/(1+0.1*y*beta))*log(1+B/(A*y*beta+Q(y*beta))) + (beta/(1.0+D))*(C+18.52*Q(beta)*D/B)*log(1+B/(A*y+10.83*Q(y*beta)+pow(y*beta,8/3.0)));

	using PhysConstants::alpha_f;

	double nu_ee = (M_e*Q(c)/hbar)*(6*pow(alpha_f,1.5)/pow(pi,2.5))*x_r*y*sqrt(beta)*I;
                                  
    //Extension to the partially degenerate case from Cassisi, S., Alexander Y. Potekhin, A. Pietrinferni, M. Catelan, and M. Salaris. 2007. “Updated Electron-Conduction Opacities: The Impact on Low-Mass Stellar Models.” The Astrophysical Journal 661 (2): 1094–1104. doi:10.1086/516819.
                                  
    double t = 25.0*T/T_F;
    double b = 135.0/sqrt(32.0*pow(pi,7.0));
                                  
    nu_ee *= (1.0+Q(t))/(1.0+t+b*Q(t)*sqrt(T/T_F));

    //Contribution due to the ei collisions according with Potekhin, A. Y., Baiko, D. A., Haensel, P., & Yakovlev, D. G. (1999). Transport properties of degenerate electrons in neutron star envelopes and white dwarf cores. Astronomy and Astrophysics, 346, 345–353. Retrieved from http://arxiv.org/abs/astro-ph/9903127

    double u_2 = 13.0;
    double u_1 = 2.8; // frequency moment of the bcc Coulomb lattice

    double Zeff = Ze[0]*X[0] + Ze[1]*X[1] + Ze[2]*X[2]  + Ze[3]*X[3] + Ze[4]*X[4] + Ze[5]*X[5]+ Ze[6]*X[6]; 
    double ai = pow(4*pi*n_i/3.0,-1.0/3.0);
	double Gamma = Q(Zeff*e_c)/(ai*k_B*T);

	beta = pi*alpha_f*Zeff*v_F/c;

	double qD = sqrt(3*Gamma)/ai;

	double k_TF2 = (alpha_f/pi)*(sqrt(1+Q(x_r))/x_r)*Q(2*k_F);

	double qi2 = Q(qD)*(1+0.06*Gamma)*exp(-sqrt(Gamma));
	double w = u_2*Q(2*k_F/qD)*(1+beta/3.0);
	double qs2 = (qi2 + k_TF2)*exp(-beta);

	double s = qs2/Q(2*k_F);

	double Lambda1, Lambda2;

	if (s < 0.005) {
		Lambda1 = 0.5*(E1(w)+log(w)+0.5772);
		Lambda2 = 0.5*(exp(-w)-1+w)/w;
	}

	else{

		if (1/w < 0.005){

			Lambda1 = 0.5*(log((1+s)/s)-(1/(s+1)));
			Lambda2 = (2*s+1)/(2*s+2) - s*log((1+s)/s);
		}

		else if (w < 0.005){

			Lambda1 = w*((2*s+1)/(2*s+2) - s*log((1+s)/s));
			Lambda2 = w*((1-3*s-6*Q(s))/(4*s+4) + 1.5*log((1+s)/s));
		}

		else {
			Lambda1 = 0.5*(log((1+s)/s)+(s/(s+1))*(1-exp(-w))-(1+s*w)*exp(s*w)*(E1(s*w)-E1(s*w+w)));
			Lambda2 = 0.5*( (exp(-w)-1+w)/w - (Q(s)/(s+1))*(1-exp(-w)) - 2*s*log((1+s)/s) +s*(2+s*w)*exp(s*w)*(E1(s*w)-E1(s*w+w)));
		}
	}

	using PhysConstants::amu;

	double Aeff = X[0]*Ai[0] + X[1]*Ai[1] + X[2]*Ai[2] + X[3]*Ai[3] + X[4]*Ai[4] + X[5]*Ai[5]+ X[6]*Ai[6]; 

	double omega_pi = sqrt(4*pi*Q(Zeff*e_c)*n_i/(Aeff*amu));

	double eta0 = 0.19/pow(Zeff,1/6.0);
	double eta = (k_B*T)/(hbar*omega_pi);
	double Gk = (eta/sqrt(Q(eta)+Q(eta0)))*(1+0.122*Q(beta)) + 0.0105*(1-1/Zeff)*(1+C(v_F/c)*beta)*eta/pow(Q(eta)+0.0081,1.5);

	double Lambda = (Lambda1-Q(v_F/c)*Lambda2)*Gk;

	double nu_ei = (4*pi*Q(Zeff*Q(e_c))*n_i/(Q(p_F)*v_F))*Lambda;

    //Calculating the total conductivity
    double nu = nu_ei + nu_ee;
                                  
    double a = Q(pi)/3.0; //strongly degenerate electron
                                  
    double kappa = a*n_e*Q(k_B)*T/(m_estar*nu);

	return log10(kappa);
}
