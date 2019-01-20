#ifndef _THERMONUCLEAR_H_
#define _THERMONUCLEAR_H_

#include <iostream>
#include <vector>
#include <string>

struct Nuclear_Network {
public:
	enum isotopes {he4, c12, o16, ne20, mg24, si28, ni56};

	enum reactions {rcag, roga, r3a, rg3a, r1212, r1216, r1616, roag, 
					rnega, rneag, rmgga, rmgag, rsiga, rcaag, rtiga,
					rsi2ni, rni2si};

	std::vector<double> rate;

	void Calculate_Rates(const double &, const double &, const std::vector<double> &);

	void TripleAlpha_Rate(); // triple alpha to 12C reaction
	void CarbonFusion_Rate(); // 12C + 12C reaction
	void CarbonOxygen_Rate(); // 12C + 16O reaction
	void OxygenFusion_Rate(); // 16O + 16O reaction
	void CarbonAlpha_Rate(); // 12C(a,g)16O and inverse
	void OxygenAlpha_Rate(); // 16O(a,g)20Ne and inverse
	void NeonAlpha_Rate(); // 20Ne(a,g)24Mg and inverse
	void MagnesiumAlpha_Rate(); // 24Mg(a,g)28Si and inverse
	void CalciumAlpha_Rate(); // 40Ca(a,g)44Ti and inverse
	void Silicon2Nickel_Rate(const std::vector<double> &); // 28Si to 56Ni and inverse

	void NuclearStatisticalEquilibrium(const double &rhop, const double &Tp, std::vector<double> &Y);

	void Electron_Screening(const std::vector<double> &Y);

	void Calculate_Derivatives(const double &, const double &,
		const std::vector<double> &, std::vector<double> &);
    void Calculate_Jacobian(const double &, const double &, const std::vector<double> &, std::vector< std::vector<double> > &);
    double EnergyGenerated(std::vector<double> &);
    double EnergyGeneratedIntegrated(const std::vector<double> &Y, const std::vector<double> &Yi, const double &dt);

private:
	double rho, T;
	double TT9, T9r, T9, T912, T913, T923, T943, T953, T932, T92, 
			T93, T972, T9r32, T9i, T9i13, T9i23, T9i32, T9i12, T9ri;
    double rsi2nida, rsi2nidsi, rni2sida;

};

#endif

