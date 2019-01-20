#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <sstream>

#include <omp.h>

#include "star.hpp"
#include "stellar_structure_equations.hpp"

#include "util/tridag.hpp"
#include "util/nr.hpp"
#include "util/roots.hpp"
#include "util/odeint.hpp"
#include "util/mins_ndim.hpp"

using namespace std;


std::string Star::ConfigFileName = "star.cfg";

std::string Star::OutputFilePath = "output/";

std::string Star::InputFilePath = "input/";

std::string Star::OutputFileName = "output";

std::string Star::EndFileName = ".dat";

std::string Star::InitialConditionFileName = "white_dwarf";

//=====================================================================
// The Shells struct
//=====================================================================

void Shells::Set_Mass(const double &Mass){
    M = Mass;
}

void Shells::Set_Temperature(const double &Temperature){
    T = Temperature;
}

void Shells::Set_Abundance(const std::vector<double> &Xin){
	X[0] = Xin[0];
    X[1] = Xin[1];
    X[2] = Xin[2];
    X[3] = Xin[3];
    X[4] = Xin[4];
    X[5] = Xin[5];
    X[6] = Xin[6];
}

void Shells::Set_Radius(const double &Radius, const double &dRadius){

    Rout = Radius;
    Rin = Radius - dRadius;
}

void Shells::Set_RadialVelocity(const double &RadialVelocity, const double &dRadialVelocity){

    Rdot_out = RadialVelocity;
    Rdot_in = RadialVelocity - dRadialVelocity;
}

void Shells::Calculate_Volume(){
    using PhysConstants::four_pi_o3;

    Vol = four_pi_o3 * (C(Rout) - C(Rin));

    if (Vol < 0){
        std::cerr << "ERROR: Shells::Calculate_Volume: Negative Volume!" << std::endl;
        exit(0);
    }

    using PhysConstants::four_pi;

    Voldot = four_pi*(Q(Rout) * Rdot_out - Q(Rin) * Rdot_in);
}

//=====================================================================
// Density within the Shell \[Rho] = M[i] / (4Pi R[i]^3/ 3)
//=====================================================================
void Shells::Calculate_Density(){
    Rho = M / Vol;

    Rhodot = -Rho * Voldot / Vol;
}
    
//=====================================================================
// Calculate the Artificial Viscosity within the Shells
//=====================================================================
void Shells::Calculate_ArtificialViscosity(){
    
    double dR = (Rout-Rin);

    const double alpha = 0.0, beta = 1.2;

    if (Rhodot > 0.0 && Rin > 0){

        double mu = dR*Rhodot/Rho;

        Q = Rho*(alpha*cs*mu + beta*Q(mu)); // the bulk viscosity as in Monaghan 1983
    }
    else Q = 0.;

}

//=====================================================================
// Calculate the material Opacity within the Shells
//=====================================================================
void Shells::Calculate_Opacity(){
    
    using PhysConstants::a;
    using PhysConstants::c;
    
    double rad = 0.2 + 1.e23*Rho*pow(T,-3.5) ;
    
    double sigmacond, kappacond;
    
    sigmacond = pow(10,cond.Evaluate_ThermalConductivity(log10(Rho),log10(T),X));
    kappacond = 4*a*c*C(T)/(3*Rho*sigmacond);

    kappa = 1./(1./kappacond+1./rad);
}

void Shells::Calculate_Ksi(){

    ksi = Ksi(Rout,Rin);

    ksidot = Ksidot(Rout,Rin,Rdot_out,Rdot_in);
}

//=====================================================================
// The Star struct
//=====================================================================

void Star::SetMass(const double &m){
    
    Mass = m;
    Distribute_Mass();
}

void Star::SetTemperature(const double &T){
    
    Temperature = T;
    Distribute_Temperature();
}

void Star::Distribute_Mass(){

	#pragma omp parallel for
    for (std::vector<int>::size_type j = 0; j < shell.size(); j++)
        shell[j].Set_Mass(Mass/Nshells);
    
}

void Star::Distribute_Radius(){
    
    using PhysConstants::four_pi_o3;
    
    shell[0].Rout = pow(shell[0].M/(four_pi_o3*Density),1/3.);
    shell[0].Rin = 0.;

    for (std::vector<int>::size_type j = 1; j < shell.size(); j++)
    {
        shell[j].Rin = shell[j-1].Rout;
        shell[j].Rout = pow(C(shell[j].Rin) + shell[j].M/(four_pi_o3*Density),1/3.);
    }   
}

void Star::Distribute_Temperature(){
    
 #pragma omp parallel for
    for (std::vector<int>::size_type j = 0; j < shell.size(); j++)
    {
        shell[j].Set_Temperature(Temperature);
    } 
}

void Star::Distribute_Abundance(){
 #pragma omp parallel for
    for (std::vector<int>::size_type j = 0; j < shell.size(); j++)
    {
        shell[j].Set_Abundance(X);
    } 
}

void Star::Initial_RadialVelocity(){
#pragma omp parallel for
    for (std::vector<int>::size_type j = 0; j < shell.size(); j++)
    {
        shell[j].Set_RadialVelocity(0.,0.);
    } 
}

void Star::Update_Radius(){

    shell[0].Rin = 0.;

    #pragma omp parallel for
    for (std::vector<int>::size_type j = 1; j < shell.size(); j++)
        shell[j].Rin = shell[j-1].Rout;
}

void Star::Update_Velocity(){

    shell[0].Rdot_in = 0.;

    #pragma omp parallel for
    for (std::vector<int>::size_type j = 1; j < shell.size(); j++){
        shell[j].Rdot_in = shell[j-1].Rdot_out ;
    }
}

void Star::Copy_InitialAbundance(){

   #pragma omp parallel for
    for (std::vector<int>::size_type j = 0; j < shell.size(); j++){
        shell[j].Xi = shell[j].X;
    }
}

void Star::Calculate_Luminosity(){

    using PhysConstants::four_pi;
    using PhysConstants::a;
    using PhysConstants::c;

    shell[0].L = four_pi*Q(shell[0].Rout)*(4*a*c/3.)*(1/(shell[0].kappa*shell[0].Rho))*pow(shell[0].T,3)*(shell[1].T-shell[0].T)/(shell[0].Rout-0.);
    
    #pragma omp parallel for
    for (std::vector<int>::size_type j = 1; j < Nshells-1; j++)
    {     
        shell[j].L = four_pi*Q(shell[j].Rout)*(4*a*c/3.)*(1/(shell[j].kappa*shell[j].Rho))*pow(shell[j].T,3)*(shell[j+1].T-shell[j].T)/(shell[j].Rout-shell[j].Rin) - four_pi*Q(shell[j-1].Rout)*(4*a*c/3.)*(1/(shell[j-1].kappa*shell[j-1].Rho))*pow(shell[j-1].T,3)*(shell[j].T-shell[j-1].T)/(shell[j-1].Rout-shell[j-1].Rin);
    }

    shell[Nshells-1].L = four_pi*Q(shell[Nshells-1].Rout)*(4*a*c/3.)*(1/(shell[Nshells-1].kappa*shell[Nshells-1].Rho))*pow(shell[Nshells-1].T,3)*(0.-shell[Nshells-1].T)/(shell[Nshells-1].Rout-shell[Nshells-1].Rin) - four_pi*Q(shell[Nshells-1].Rin)*(4*a*c/3.)*(1/(shell[Nshells-2].kappa*shell[Nshells-2].Rho))*pow(shell[Nshells-2].T,3)*(shell[Nshells-1].T-shell[Nshells-2].T)/(shell[Nshells-2].Rout-shell[Nshells-2].Rin);
    
    Luminosity = - four_pi*Q(shell[Nshells-1].Rout)*(4*a*c/3.)*(1/(shell[Nshells-1].kappa*shell[Nshells-1].Rho))*pow(shell[Nshells-1].T,3)*(0.-shell[Nshells-1].T)/(shell[Nshells-1].Rout-shell[Nshells-1].Rin) ;
}


void Star::Update_Interior(){

    Update_Radius();
    Update_Velocity();

    for (std::vector<int>::size_type j = 0; j < shell.size(); j++){
        shell[j].Calculate_Ksi();
        
        shell[j].Calculate_Volume();
        shell[j].Calculate_Density();

        eos.Evaluate(shell[j].Rho,shell[j].T,shell[j].X);

        shell[j].P = eos.P;
        shell[j].u = eos.u;
    
        shell[j].cs = eos.cs;
    
        shell[j].Cv = shell[j].M * fabs(eos.cv);
        shell[j].Cp = shell[j].M * fabs(eos.cp);
        shell[j].dUdrho = shell[j].M * eos.dudrho;

        if (art_viscosity) shell[j].Calculate_ArtificialViscosity();
    
        if (luminosity) shell[j].Calculate_Opacity();
    }
        
    if (luminosity) Calculate_Luminosity();
}

void Star::ReadConfigFile(){
    
    string ReadFileName = InputFilePath + ConfigFileName;
    
    ifstream InFile(ReadFileName.c_str());
    
    std::cerr << "Reading configuration file: " << ConfigFileName << " ...";
    
    if (!InFile.is_open()) {
        std::cout << "The configuration file could not be opened (wrong filename?)" << std::endl;
        exit(1);
    }
    
    else{
        std::string line;
        while(getline(InFile, line)){
            // read all lines until end of configuration file is reached
            if ( line[0]=='#' || line[0]==' ' || line[0]=='\n' ) {} /* Ignore lines beginning with empty space, #
                                                                     or newline character. */
            else{
                if (line.find(" ")== std::string::npos) {
                    // if there are no empty spaces, we are done
                }
                else {
                    line = line.substr(0,line.find_first_of(" "));  // remove everything after first empty space
                }
                
                ProcessLine(&line); // use helper function to sort into respective variables
            }
        }

        InFile.close();
    }
    
    InFile.clear();
    
    std::cerr << "done."<< std::endl;
    
}

void Star::ProcessLine(std::string *line){
    
    std::string pre = line->substr(0,line->find("="));   // part of the string in front of "="
    std::string post = line->substr(line->find("=")+1);  // part of the string after "="
    
    if (pre == "Nshells") Nshells = atof(post.c_str()); 	// star initial temperature
    
    if (pre == "eos") eos.id_eos = atoi(post.c_str());
    
    if (pre == "Gravity") {
        if (post == "yes") gravity = true;
        else if (post == "no") gravity = false;
    }

    if (pre == "Reactions") {
        if (post == "yes") reactions = true;
        else if (post == "no") reactions = false;
    }

    if (pre == "Art_Viscosity") {
        if (post == "yes") art_viscosity = true;
        else if (post == "no") art_viscosity = false;
    }

    if (pre == "Luminosity") {
        if (post == "yes") luminosity = true;
        else if (post == "no") luminosity = false;
    }

    if (pre == "OutputFileName") OutputFileName = post;
    
    if (pre == "OutputFilePath") OutputFilePath = post;
    
    if (pre == "InitialConditionFileName") InitialConditionFileName = post;
    
}


void Star::Read_InitialCondition(){
    
    string ReadFileName = InputFilePath + InitialConditionFileName;
    
    ifstream InFile(ReadFileName.c_str());
    
    std::cerr << "Reading the initial condition file: " << InitialConditionFileName << "...";
    
    if (!InFile) {
        std::cout << "The initial condition file could not be opened (wrong filename?)" << std::endl;
        exit(1);
    }
    
    std::string line;
    unsigned int j = 0;

    if (InFile.is_open()){
        getline(InFile, line);  // read all lines until end of configuration file is reached
        getline(InFile, line);
    }

    Mass = 0.;
    while(!InFile.eof()){
        InFile >> shell[j].M >> shell[j].Rout >> shell[j].Rdot_out >> shell[j].T >> shell[j].X[0] >> shell[j].X[1] >> shell[j].X[2] >> shell[j].X[3] >> shell[j].X[4] >> shell[j].X[5] >> shell[j].X[6];
            j++;
        Mass += shell[j].M;
    }

    InFile.close();
    
    std::cerr << "done."<< std::endl;
    
}

//=====================================================================
// Função que calcula a energia gravitacional Vg total da estrela
// para uma dada configuração {R_i}
//=====================================================================

double Star::Gravitational_Energy(){

    using PhysConstants::G;
    const double three_teenths = 0.3;

    std::vector<double> sum_Mass(Nshells-1); //Vetor para a soma das massas das camadas internas

    sum_Mass[0] = shell[0].M;
    for (std::vector<int>::size_type j = 1; j < Nshells-1; j++)
        sum_Mass[j] = sum_Mass[j-1] + shell[j].M;

    shell[0].Vg = -G*three_teenths * (shell[0].M / shell[0].Rout)
            * (f(shell[0].ksi) * shell[0].M ); // For first shell
    
    #pragma omp parallel for
    for (std::vector<int>::size_type j = 1; j < Nshells; j++){ 
        shell[j].Vg = -G*three_teenths * (shell[j].M / shell[j].Rout)
            * (f(shell[j].ksi) * shell[j].M + g(shell[j].ksi) * (sum_Mass[j-1]));
    }

    double sum = 0.0;
    for (std::vector<int>::size_type j = 0; j < Nshells; j++)
        sum += shell[j].Vg;

    sum_Mass.clear();
    
    return (sum);
}

//=====================================================================
// Função que calcula a energia interna U_int total da estrela
// para uma dada configuração {R_i}
//=====================================================================

double Star::Internal_Energy() {

    #pragma omp parallel for
    for (std::vector<int>::size_type j = 0; j < Nshells; j++){ 
        shell[j].U = shell[j].u*shell[j].M;
    }

    double sum = 0.0;
    for (std::vector<int>::size_type j = 0; j < Nshells; j++) {
        sum += shell[j].U ;
    }

    return sum;
}

double Star::Nuclear_Energy() {

    double sum = 0.;
    for (std::vector<int>::size_type j = 0; j < Nshells; j++) {
        sum += shell[j].M*shell[j].network.EnergyGeneratedIntegrated(shell[j].X,shell[j].Xi,1.);
    }

    return sum;
}

//=====================================================================
// Função que calcula a energia cinética K total da estrela
// para uma dada configuração {R_i,Rdot_i}
//=====================================================================

double Star::Kinetic_Energy() {

    const double half = 0.5;

    double K;

    std::vector<double> Td(Nshells), Ts(Nshells);

    T_matrix(Td,Ts);

    if (Nshells==1) K = half*shell[0].Rdot_out*Td[0]*shell[0].Rdot_out;
    
    else
    {
        double sum = 0.0;
        
        for (std::vector<int>::size_type j = 0; j < Nshells-1; j++){
            sum += Td[j]*Q(shell[j].Rdot_out) + 2*Ts[j]*shell[j].Rdot_out*shell[j+1].Rdot_out;
        }
        sum += Td[Nshells-1]*Q(shell[Nshells-1].Rdot_out);

        K = half*sum;
    }

    Td.clear(); Ts.clear();

    return K;
}

//=====================================================================
// Matriz tridiagonal de energia cinética n x n (T) a partir dos
// vetores a, b das diagonais
//=====================================================================

void Star::T_matrix(vector<double>& d, vector<double>& s)
{
    if (Nshells==1) d[0] = T_22(shell[0].ksi, shell[0].M);

    else{
    // Create the a, b and c coefficients in the vectors
       #pragma omp parallel for
        for (std::vector<int>::size_type j = 0; j < Nshells-1; j++)
        {
          d[j] = T_22(shell[j].ksi, shell[j].M) + T_11(shell[j+1].ksi, shell[j+1].M);
          s[j] = T_12(shell[j+1].ksi, shell[j+1].M);
        }

        d[Nshells-1] = T_22(shell[Nshells-1].ksi, shell[Nshells-1].M);
        s[Nshells-1] = 0.; // por definição de matriz tridiagonal simétrica
    }
}

//=====================================================================
// Matriz tridiagonal de energia cinética n x n (Q) a partir dos
// vetores a, b e c das diagonais
//=====================================================================

void Star::Q_matrix(vector<double>& d, vector<double>& s)
{
    if (Nshells==1) d[0] = Q_22(shell[0].ksi, shell[0].ksidot, shell[0].M);

    else{ 
    // Create the a, b and c coefficients in the vectors
      #pragma omp parallel for
        for (std::vector<int>::size_type j = 0; j < shell.size()-1; j++)
        {
          d[j] = Q_22(shell[j].ksi, shell[j].ksidot, shell[j].M) + Q_11(shell[j+1].ksi, shell[j+1].ksidot, shell[j+1].M);
          s[j] = Q_12(shell[j+1].ksi, shell[j+1].ksidot, shell[j+1].M);
        }

      d[Nshells-1] = Q_22(shell[Nshells-1].ksi, shell[Nshells-1].ksidot, shell[Nshells-1].M);
      s[Nshells-1] = 0.; // por definição de matriz tridiagonal simétrica
    }
}


//=====================================================================
// Vetor de força 1 x n (F)
//=====================================================================

void Star::F_vector(std::vector<double>& f)
{
    using PhysConstants::four_pi;

    if (gravity){
        using PhysConstants::G;

        if (Nshells==1) f[0] = four_pi*Q(shell[0].Rout)*( shell[0].P + shell[0].Q ) - (G/Q(shell[0].Rout))*( f_1(shell[0].ksi)*Q(shell[0].M) );

        else{

            std::vector<double> sumM(Nshells-1); //Vetor para a soma das massas das camadas internas

            sumM[0] = shell[0].M;
            for (std::vector<int>::size_type j = 1; j < shell.size()-1; j++)
                sumM[j] = sumM[j-1] + shell[j].M;

            f[0] = four_pi*Q(shell[0].Rout)*( (shell[0].P-shell[1].P) + (shell[0].Q-shell[1].Q) ) - (G/Q(shell[0].Rout))*( f_1(shell[0].ksi)*Q(shell[0].M) + f_2(shell[1].ksi)*Q(shell[1].M) + g_2(shell[1].ksi)*shell[1].M*sumM[0] );

            // A seguir, calculamos as variáveis necessárias para o elemento F[i] dentro do loop
         #pragma omp parallel for
            for (std::vector<int>::size_type j = 1; j < shell.size()-1; j++)
            {
                f[j] = four_pi*Q(shell[j].Rout)*( (shell[j].P-shell[j+1].P) + (shell[j].Q-shell[j+1].Q) ) - (G/Q(shell[j].Rout))*( f_1(shell[j].ksi)*Q(shell[j].M) + f_2(shell[j+1].ksi)*Q(shell[j+1].M) + g_1(shell[j].ksi)*shell[j].M*sumM[j-1]  + g_2(shell[j+1].ksi)*shell[j+1].M*sumM[j] );
            }
        
            // A seguir, calculamos as variáveis necessárias para o elemento F[N-1]
            f[Nshells-1] = four_pi*Q(shell[Nshells-1].Rout)*(shell[Nshells-1].P + shell[Nshells-1].Q) - (G/Q(shell[Nshells-1].Rout))*( f_1(shell[Nshells-1].ksi)*Q(shell[Nshells-1].M) + g_1(shell[Nshells-1].ksi)*shell[Nshells-1].M*sumM[Nshells-2] );
        
            sumM.clear();
        }

    }

    else {

        if (Nshells==1) f[0] = four_pi*Q(shell[0].Rout)*( shell[0].P + shell[0].Q );

        else{

            f[0] = four_pi*Q(shell[0].Rout)*( (shell[0].P-shell[1].P) + (shell[0].Q-shell[1].Q) );

            // A seguir, calculamos as variáveis necessárias para o elemento F[i] dentro do loop
         #pragma omp parallel for
            for (std::vector<int>::size_type j = 1; j < shell.size()-1; j++)
            {
                f[j] = four_pi*Q(shell[j].Rout)*( (shell[j].P-shell[j+1].P) + (shell[j].Q-shell[j+1].Q) );
            }
        
            // A seguir, calculamos as variáveis necessárias para o elemento F[N-1]
            f[Nshells-1] = four_pi*Q(shell[Nshells-1].Rout)*(shell[Nshells-1].P + shell[Nshells-1].Q);
        }

    }
}

//=====================================================================
// Retorna as derivadas dR/dt (Velocidade) para a dinâmica
//=====================================================================

void Star::dRdt()
{  
  #pragma omp parallel for
    for (std::vector<int>::size_type j = 0; j < shell.size(); j++)
        shell[j].dRdt = shell[j].Rdot_out; 
}

//=====================================================================
// Retorna as derivadas dRdot/dt (Aceleração) para a dinâmica
//=====================================================================

void Star::dRdotdt()
{    
    std::vector<double> Td(Nshells), Ts(Nshells);    
    std::vector<double> Qd(Nshells), Qs(Nshells);
    std::vector<double> f(Nshells), aux(Nshells), fout(Nshells);

    T_matrix(Td, Ts);
    Q_matrix(Qd, Qs);
    F_vector(f);

    if (Nshells==1) fout[0] = shell[0].Rdot_out*Qd[0] + f[0];
    
    else
    {
        aux[0] = shell[0].Rdot_out*Qd[0] + shell[1].Rdot_out*Qs[0] + f[0];
        
       #pragma omp parallel for
        for (std::vector<int>::size_type j = 1; j < shell.size()-1; j++)
            aux[j] = shell[j-1].Rdot_out*Qs[j-1] + shell[j].Rdot_out*Qd[j] + shell[j+1].Rdot_out*Qs[j] + f[j];

        aux[Nshells-1] = shell[Nshells-1].Rdot_out*Qd[Nshells-1] + shell[Nshells-2].Rdot_out*Qs[Nshells-2] + f[Nshells-1];

        tridag_symm(Ts,Td,aux,fout);
    }

    #pragma omp parallel for
    for (std::vector<int>::size_type j = 0; j < shell.size(); j++)
        shell[j].dRdotdt = fout[j];

    Td.clear(); Ts.clear(); Qd.clear(); Qs.clear(); f.clear(); aux.clear(); fout.clear();
}

//=====================================================================
// Retorna as derivadas dT/dt para a dinâmica
//=====================================================================
void Star::ThermalConductionEquation()
{
  #pragma omp parallel for
    for (std::vector<int>::size_type j = 0; j < shell.size(); j++){
        if (tau > 0.0) shell[j].dTdt = (-shell[j].dUdrho*shell[j].Rhodot - shell[j].P*shell[j].Voldot +shell[j].Qd + shell[j].L ) /(shell[j].Cv);
        
        else shell[j].dTdt = (-shell[j].dUdrho*shell[j].Rhodot - shell[j].P*shell[j].Voldot -shell[j].Q*shell[j].Voldot + shell[j].L + shell[j].Eps*shell[j].M + shell[j].dUdt) /(shell[j].Cv);
    }
}

//=====================================================================
// Retorna as derivadas dQd/dt (calor retardado) para a dinâmica
//=====================================================================
void Star::dQddt()
{  
  #pragma omp parallel for
    for (std::vector<int>::size_type j = 0; j < shell.size(); j++){
        shell[j].dQddt = (-shell[j].Q*shell[j].Voldot + shell[j].Eps*shell[j].M -shell[j].Qd)/tau;
    }
}

//=====================================================================
// Funções para cálculo do equilíbrio hidrostático da estrela
//=====================================================================

struct Equilibrium {
    Star *star_ptr;
    
    void initial(std::vector<double> &x){
        
        #pragma omp parallel for
        for (std::vector<int>::size_type j = 0; j < star_ptr->Nshells; j++){
            x[j]= sqrt(star_ptr->shell[j].Rout - star_ptr->shell[j].Rin);
        }
    }
    
    double operator()(std::vector<double> &x){

        using PhysConstants::four_pi_o3;
        const double onethird = 1.0/3.0;

        star_ptr->shell[0].Rout = Q(x[0]);
        star_ptr->shell[0].Rin = 0.;
        for (std::vector<int>::size_type j = 1; j < star_ptr->Nshells; j++) {
            star_ptr->shell[j].Rin = star_ptr->shell[j-1].Rout;
            star_ptr->shell[j].Rout = star_ptr->shell[j].Rin + Q(x[j]);
        }
        
        star_ptr->Update_Interior();

        std::cerr << "Energy = " << (star_ptr->Gravitational_Energy() + star_ptr->Internal_Energy()) << std::endl;
        
        return (star_ptr->Gravitational_Energy() + star_ptr->Internal_Energy());
    }
    
};


void Star::HydrostaticEquilibrium(const double &rtol){
    
    std::vector<double> x(Nshells);
    
    Equilibrium energy;
    energy.star_ptr = this;
    energy.initial(x);
    
    Powell<Equilibrium> powell(energy,rtol);
    
    x = powell.minimize(x);
    
    x.clear();
}

struct Network: public Nuclear_Network {
    double Rho, T;
    
    void operator() (const double &t, const std::vector<double> &y, std::vector<double> &dydt)
    {
        Calculate_Derivatives(Rho,T,y,dydt);
    }
    
    void jacobian(const double &t, const std::vector<double> &y, std::vector<double> &dydt, std::vector< std::vector<double> > &dfdy)
    {
        Calculate_Jacobian(Rho,T,y,dfdy);
    }
};


void Star::NuclearReactions(double &dtnet, const double &rtol, const double &tol){
    
    double hmin = 0.0, h0 = 1.0e-9;

    std::vector<int>::size_type j;

    //#pragma omp parallel for private(j,hmin,h0,rtol,tol,dtnet)
    for (j = 0; j < shell.size(); j++){
        
        Network network;
      
        std::vector<double> Xi = shell[j].X;
      
        network.Rho = shell[j].Rho;
        network.T = shell[j].T;
        
        double dtNSE = pow(shell[j].Rho,0.2)*exp(1.797e11/shell[j].T-40.5);

        double h = h0;
        double t = 0.0;

//        if (dtNSE <= dtnet){
//            network.NuclearStatisticalEquilibrium(shell[j].Rho,shell[j].T,shell[j].X);
//        }
//
//        else{
            Odeint<StepperSie<Network> > odenet(shell[j].X,t,dtnet,tol,rtol,h,hmin,network);
                
            odenet.integrate();
//        }

        shell[j].Eps = network.EnergyGeneratedIntegrated(shell[j].X,Xi,dtnet);
        
        Xi.clear();
    }
  
}

void Star::Single_Step(const double &dt, const double &tol, const double &rtol){
    
    Update_Interior();
    dRdt();
    dRdotdt();
    ThermalConductionEquation();
    if (tau > 0.0) dQddt();
}


double Star::TimeStep(unsigned int id_dynamics){

    double dtf, dtR, dtacc, dtT, dtQd, dtz, dtcv, dtdyn, dtexpl;

    if (id_dynamics == 1){
        double dR = (shell[0].Rout-shell[0].Rin);
        dtacc = sqrt(dR/fabs(shell[0].dRdotdt));
        dtR = shell[0].Rout/fabs(shell[0].Rdot_out);
        dtcv = dR/(shell[0].cs);
        dtdyn = 446/sqrt(shell[0].Rho);

        for (std::vector<int>::size_type j = 1; j < shell.size(); j++){
            dR = (shell[j].Rout-shell[j].Rin);
            dtR = MIN(dtR,shell[j].Rout/fabs(shell[j].Rdot_out));
            dtacc = sqrt(dR/fabs(shell[j].dRdotdt));
            dtdyn = 446/sqrt(shell[j].Rho);
            dtcv = dR/(shell[j].cs+dR*fabs(shell[j].Rhodot/shell[j].Rho));
        }

        std::cerr << "dtR = " << dtR << std::endl;
        std::cerr << "dtacc = " << dtacc << std::endl;  
        std::cerr << "dtdyn = " << dtdyn << std::endl;
        std::cerr << "dtcv = " << dtcv << std::endl;

        return 0.25*MIN(dtcv,MIN(dtdyn,MIN(dtR,dtacc)));
    }

    if (id_dynamics == 2){
        dtT = fabs(shell[0].T/fabs(shell[0].dTdt));
        dtexpl = fabs((shell[0].Cp*shell[0].T)/(shell[0].Eps*shell[0].M));

        for (std::vector<int>::size_type j = 1; j < shell.size(); j++){
            dtT = MIN(dtT,fabs(shell[j].T/shell[j].dTdt));
            dtexpl = MIN(dtexpl,fabs((shell[j].Cp*shell[j].T)/(shell[j].Eps*shell[j].M)));
        }

        return 0.25*MIN(dtT,dtexpl);
    }

    if (id_dynamics > 2 ){
        double dR = (shell[0].Rout-shell[0].Rin);
        dtdyn = 446/sqrt(shell[0].Rho);
        dtR = shell[0].Rout/fabs(shell[0].Rdot_out);
        dtacc = sqrt(dR/fabs(shell[0].dRdotdt));
        dtcv = dR/(shell[0].cs);
        dtT = shell[0].T/fabs(shell[0].dTdt);
        dtexpl = fabs((shell[0].Cp*shell[0].T)/(shell[0].Eps*shell[0].M));
        double dtNSE = pow(shell[0].Rho,0.2)*exp(1.797e11/shell[0].T-40.5);
        
        for (std::vector<int>::size_type j = 1; j < shell.size(); j++){
            dR = (shell[j].Rout-shell[j].Rin);
            dtdyn = MIN(dtdyn,446/sqrt(shell[j].Rho));
            dtcv = MIN(dtcv,dR/(shell[j].cs));
            dtexpl = MIN(dtexpl,fabs((shell[j].Cp*shell[j].T)/(shell[j].Eps*shell[j].M)));
            dtR = MIN(dtR,shell[j].Rout/fabs(shell[j].Rdot_out));
            dtacc = sqrt(dR/fabs(shell[j].dRdotdt));
            dtT = MIN(dtT,shell[j].T/fabs(shell[j].dTdt));
            
            dtNSE = MIN(dtNSE,pow(shell[j].Rho,0.2)*exp(1.797e11/shell[j].T-40.5));
        }

        std::cerr << "dtR = " << dtR << std::endl;
        std::cerr << "dtT = " << dtT << std::endl;
        std::cerr << "dtacc = " << dtacc << std::endl;
        std::cerr << "dtdyn = " << dtdyn << std::endl;
        std::cerr << "dtcv = " << dtcv << std::endl;
        std::cerr << "dtexpl = " << dtexpl << std::endl;
        std::cerr << "dtNSE = " << dtNSE << std::endl;
        
        return 0.3*MIN(dtNSE,MIN(dtT,MIN(dtR,MIN(dtacc,MIN(dtexpl,MIN(dtdyn,dtcv))))));
    }
    
}

struct Dynamics{
    Star *starptr;
    unsigned int id;
    double ti;
    
    void operator() (const double &t, const std::vector<double> &y, std::vector<double> &dydt)
    {
        unsigned int n = starptr->shell.size();
        
        if (id == 1){
            
#pragma omp parallel for
            for (std::vector<int>::size_type j = 0; j < starptr->shell.size(); j++){
                starptr->shell[j].Rout = y[j];
                starptr->shell[j].Rdot_out = y[j+n];
            }
            
            starptr->Update_Interior();
            starptr->dRdt();
            starptr->dRdotdt();
            
            #pragma omp parallel for
            for (std::vector<int>::size_type j = 0; j < starptr->shell.size(); j++){
                dydt[j] = starptr->shell[j].dRdt;
                dydt[j+n] = starptr->shell[j].dRdotdt;
            }
        }
        
        if (id == 2){
            #pragma omp parallel for
            for (std::vector<int>::size_type j = 0; j < starptr->shell.size(); j++){
                starptr->shell[j].T = y[j];
            }
            
            starptr->Update_Interior();
            starptr->ThermalConductionEquation();
            
            #pragma omp parallel for
            for (std::vector<int>::size_type j = 0; j < starptr->shell.size(); j++){
                dydt[j] = starptr->shell[j].dTdt;
            }
            
        }
        
        if (id == 3){

        #pragma omp parallel for
            for (std::vector<int>::size_type j = 0; j < starptr->shell.size(); j++){
                starptr->shell[j].Rout = y[j];
                starptr->shell[j].Rdot_out = y[j+n];
                starptr->shell[j].T = y[j+2*n];
                if (starptr->tau > 0.0) starptr->shell[j].Qd = y[j+3*n];
            }
            
            starptr->Update_Interior();
            starptr->dRdt();
            starptr->dRdotdt();
            
            starptr->ThermalConductionEquation();
            if (starptr->tau > 0.0) starptr->dQddt();
            
        #pragma omp parallel for
            for (std::vector<int>::size_type j = 0; j < starptr->shell.size(); j++){
                dydt[j] = starptr->shell[j].dRdt;
                dydt[j+n] = starptr->shell[j].dRdotdt;
                dydt[j+2*n] = starptr->shell[j].dTdt;
                if (starptr->tau > 0.0) dydt[j+3*n] = starptr->shell[j].dQddt;
            }
            
        }
        
    }
        
};


void Star::Time_Integration(double &ti, double &dt, const double &rtol, const double &tol, unsigned int id_dynamics){
    
    double hmin = 0.;
    double dtnet;
    
    double tf = ti + dt;
    
    Dynamics dynamics;
    dynamics.starptr = this;
    dynamics.id = id_dynamics;
    dynamics.ti = ti;
    
    const unsigned int n = shell.size();
    
    if (id_dynamics == 1){
        
        std::vector<double> y(2*n);
        
    #pragma omp parallel for
        for (std::vector<int>::size_type j = 0; j < shell.size(); j++){
            y[j] = shell[j].Rout;
            y[j+n] = shell[j].Rdot_out;
        }
        
        Odeint<StepperFehlb<Dynamics> > ode(y,ti,tf,tol,rtol,dt_integrator,hmin,dynamics);
        
        ode.integrate();
        
        dt_integrator = ode.s.hnext; //Atualiza o tamanho do passo
        
        y.clear();
        
    }
    
    if (id_dynamics == 2){
        
        std::vector<double> y(n);
        
        #pragma omp parallel for
        for (std::vector<int>::size_type j = 0; j < shell.size(); j++){
            y[j] = shell[j].T;
        }
        
        Odeint<StepperFehlb<Dynamics> > ode(y,ti,tf,tol,rtol,dt_integrator,hmin,dynamics);
        
        ode.integrate();
        
        dt_integrator = ode.s.hnext; //Atualiza o tamanho do passo
        
        y.clear();
    }
    
    if (id_dynamics == 3){
        
        dtnet = dt;
        
        Update_Interior();
        if (reactions) NuclearReactions(dtnet,1.0e-6,1.0e-8);
        
        std::vector<double> y(3*n);
        if (tau > 0.0) y.resize(4*n);
        
        #pragma omp parallel for
        for (std::vector<int>::size_type j = 0; j < shell.size(); j++){
            y[j] = shell[j].Rout;
            y[j+n] = shell[j].Rdot_out;
            y[j+2*n] = shell[j].T;
            if (tau > 0.0) y[j+3*n] = shell[j].Qd;
        }
        
        Odeint<StepperFehlb<Dynamics> > ode(y,ti,tf,tol,rtol,dt_integrator,hmin,dynamics);
        
        ode.integrate();
        
        dt_integrator = ode.s.hnext; //Atualiza o tamanho do passo
        
        y.clear();
    }
    
    dt = MIN(10*dt_integrator,TimeStep(id_dynamics));
    
    std::cerr << "dt =" << dt << std::endl;
    
    dt_integrator = MIN(dt_integrator,dt);

    std::cerr << "dt_integrator = " << dt_integrator << std::endl;
    
    ti = tf;
}
