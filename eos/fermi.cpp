#include "fermi.hpp"
#include "../util/quadrature.hpp"
#include "../util/derule.hpp"

double FermiDirac::operator() (const double &t){
	//Integrando para quadratura trapezoidal da integral generalizada de Fermi-Dirac com transformação x = exp(t-e^{-t})
	double x;
	x = exp(t-exp(-t));
	return x*(1.0+exp(-t))*pow(x,kk)*sqrt(1.0+thetaa*0.5*x)/(exp(x-etaa)+1.0);
}

double FermiDirac::operator() (const double &x, const double &del){
	//Integrando para DE de quadratura da integral de Fermi-Dirac generalizada
	if (x < 1.0)
		return pow(del,kk)*sqrt(1.0+thetaa*0.5*x)/(exp(x-etaa)+1.0);
	else 		
		return pow(x,kk)*sqrt(1.0+thetaa*0.5*x)/(exp(x-etaa)+1.0);

}

double FermiDirac::F(const double &k, const double &eta, const double &theta)
	//Computa a integral de Fermi-Dirac generalizada F_k(eta,theta) onde k > -1 e \theta >= 0. A acurácia é aproximada-
	// mente o quadrado do parâmetro EPS. NMAX limita o número total de passos da quadratura.
{
	const double EPS = 1.e-8;
	const double TINY = std::numeric_limits<double>::epsilon();
	const int NMAX = 16;
	double a, aa, b, bb, hmax, olds, sum;
	kk=k; // Armazena os argumentos nas variáveis membros para uso nas avaliações das funções
	etaa=eta;
	thetaa=theta;

	//if (etaa < -1/thetaa*EPS) return 0.;

	if(eta <= 15.0){
		a=-4.5;		//Atribua limites para o mapeamento x=exp(t-e^{-t})
		b=5.0;
		Trapzd<FermiDirac> s(*this,a,b);
		for (unsigned int i=1;i<=NMAX;i++){
			sum=s.next();
			if (i>3){ 		//Teste para convergência
				if(fabs(sum-olds) <= EPS*fabs(olds) + TINY)
					return sum;
			}
			olds=sum; //Valor salvo para próximo teste de convergência
		}
	}
	else {
		a=0.0;			//Designe limites para a regra DE
		b=eta;
		aa=b;
		bb=eta+60.0;
		hmax=4.3;		//Grande o suficiente para controlar k negativo ou \eta grande
		DErule<FermiDirac> s(*this,a,b,hmax);
		DErule<FermiDirac> ss(*this,aa,bb,hmax);
		for (unsigned int i=1;i<=NMAX;i++){
			sum=s.next()+ss.next();
			if(i>3)
				if(fabs(sum-olds) <= EPS*fabs(olds))
					return sum;
			olds=sum;
		}
	}
	std::cerr << "FermiDirac function: No convergence in eta = " << eta << " and beta = " << theta << std::endl;
	return 0.0;
}

double FermiDiracDBeta::operator() (const double &t){
	//Integrando para quadratura trapezoidal da integral generalizada de Fermi-Dirac com transformação x = exp(t-e^{-t})
	double x;
	x = exp(t-exp(-t));
	return 0.25*x*(1.0+exp(-t))*pow(x,kk+1)/(sqrt(1.0+thetaa*0.5*x)*(exp(x-etaa)+1.0) );
}

double FermiDiracDBeta::operator() (const double &x, const double &del){
	//Integrando para DE de quadratura da integral de Fermi-Dirac generalizada
	if (x < 1.0)
		return 0.25*pow(del,kk+1)/(sqrt(1.0+thetaa*0.5*x)*(exp(x-etaa)+1.0) );
	else 		
		return 0.25*pow(x,kk+1)/(sqrt(1.0+thetaa*0.5*x)*(exp(x-etaa)+1.0) );

}

double FermiDiracDBeta::dFdbeta(const double &k, const double &eta, const double &theta)
	//Computa a integral de Fermi-Dirac generalizada F_k(eta,theta) onde k > -1 e \theta >= 0. A acurácia é aproximada-
	// mente o quadrado do parâmetro EPS. NMAX limita o número total de passos da quadratura.
{
	const double EPS = 1.e-8;
	const double TINY = std::numeric_limits<double>::epsilon();
	const int NMAX = 16;
	double a, aa, b, bb, hmax, olds, sum;
	kk=k; // Armazena os argumentos nas variáveis membros para uso nas avaliações das funções
	etaa=eta;
	thetaa=theta;

	if(eta <= 15.0){
		a=-4.5;		//Atribua limites para o mapeamento x=exp(t-e^{-t})
		b=5.0;
		Trapzd<FermiDiracDBeta> s(*this,a,b);
		for (unsigned int i=1;i<=NMAX;i++){
			sum=s.next();
			if (i>3){ 		//Teste para convergência
				if(fabs(sum-olds) <= EPS*fabs(olds) + TINY)
					return sum;
			}
			olds=sum; //Valor salvo para próximo teste de convergência
		}
	}
	else {
		a=0.0;			//Designe limites para a regra DE
		b=eta;
		aa=b;
		bb=eta+60.0;
		hmax=4.3;		//Grande o suficiente para controlar k negativo ou \eta grande
		DErule<FermiDiracDBeta> s(*this,a,b,hmax);
		DErule<FermiDiracDBeta> ss(*this,aa,bb,hmax);
		for (unsigned int i=1;i<=NMAX;i++){
			sum=s.next()+ss.next();
			if(i>3)
				if(fabs(sum-olds) <= EPS*fabs(olds))
					return sum;
			olds=sum;
		}
	}
	std::cerr << "FermiDiracDBeta function: No convergence in eta = " << eta << " and beta = " << theta << std::endl;
	return 0.0;
}



double FermiDiracDEta::operator() (const double &t){
	//Integrando para quadratura trapezoidal da integral generalizada de Fermi-Dirac com transformação x = exp(t-e^{-t})
	double x = exp(t-exp(-t));
	return x*(1.0+exp(-t))*pow(x,kk)*sqrt(1.0+thetaa*0.5*x)/SQR(exp((x-etaa)/2)+exp((etaa-x)/2));
}

double FermiDiracDEta::operator() (const double &x, const double &del){
	//Integrando para DE de quadratura da integral de Fermi-Dirac generalizada
	if (x < 1.0)
		return pow(del,kk)*sqrt(1.0+thetaa*0.5*x)/SQR(exp((x-etaa)/2)+exp((etaa-x)/2));
	else 		
		return pow(x,kk)*sqrt(1.0+thetaa*0.5*x)/SQR(exp((x-etaa)/2)+exp((etaa-x)/2));

}

double FermiDiracDEta::dFdeta(const double &k, const double &eta, const double &theta)
	//Computa a integral de Fermi-Dirac generalizada F_k(eta,theta) onde k > -1 e \theta >= 0. A acurácia é aproximada-
	// mente o quadrado do parâmetro EPS. NMAX limita o número total de passos da quadratura.
{
	const double EPS = 1.e-8;
	const double TINY = std::numeric_limits<double>::epsilon();
	const int NMAX = 16;
	double a, aa, b, bb, hmax, olds, sum;
	kk=k; // Armazena os argumentos nas variáveis membros para uso nas avaliações das funções
	etaa=eta;
	thetaa=theta;

	if(eta <= 15.0){
		a=-4.5;		//Atribua limites para o mapeamento x=exp(t-e^{-t})
		b=5.0;
		Trapzd<FermiDiracDEta> s(*this,a,b);
		for (unsigned int i=1;i<=NMAX;i++){
			sum=s.next();
			if (i>3){ 		//Teste para convergência
				if(fabs(sum-olds) <= EPS*fabs(olds) + TINY)
					return sum;
			}
			olds=sum; //Valor salvo para próximo teste de convergência
		}
	}
	else {
		a=0.0;			//Designe limites para a regra DE
		b=eta;
		aa=b;
		bb=eta+60.0;
		hmax=4.3;		//Grande o suficiente para controlar k negativo ou \eta grande
		DErule<FermiDiracDEta> s(*this,a,b,hmax);
		DErule<FermiDiracDEta> ss(*this,aa,bb,hmax);
		for (unsigned int i=1;i<=NMAX;i++){
			sum=s.next()+ss.next();
			if(i>3)
				if(fabs(sum-olds) <= EPS*fabs(olds))
					return sum;
			olds=sum;
		}
	}
	std::cerr << "FermiDiracDeta function: No convergence in eta = " << eta << " and beta = " << theta << std::endl;
	return 0.0;
}
