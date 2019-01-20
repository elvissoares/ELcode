/*user_IO

    General Description:
    ====================
        This file takes care of all the user input requests and most output
---------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "cpp_constants_CGS.hpp"
#include "star.hpp"
#include "eos.hpp"

#define tab "\t"

struct OutDynamics {

    Star *star_ptr;

    std::ofstream OutRadius, OutVelocity, OutMass, OutDensity, 
                    OutTemperature, OutPressure, OutInternalEnergy, OutHelium, OutCarbon,
                    OutOxygen, OutNeon, OutMagnesium, OutSilicon, OutNickel, OutNuclearEnergy,
                    OutEnergy, OutLuminosity, OutSoundSpeed, OutCentralInfo;

    OutDynamics(Star &star);

    void RadiusFile();
    void VelocityFile();
    void MassFile();
    void DensityFile(); 
    void TemperatureFile();
    void InternalEnergyFile();
    void PressureFile(); 
    void SoundSpeedFile(); 

    void HeliumFile();
    void CarbonFile();
    void OxygenFile();   
    void NeonFile();
    void MagnesiumFile();   
    void SiliconFile();
    void NickelFile();
    void NuclearEnergyFile();

    void EnergyFile();

    void Luminosity();

    void AllVariables(const double &t);

    void CentralInfo();

};

OutDynamics::OutDynamics(Star &star)
{
    star_ptr = &star;       
    RadiusFile();
    VelocityFile();
    MassFile();
    DensityFile();
    TemperatureFile();
    PressureFile();
    InternalEnergyFile();
    SoundSpeedFile(); 
    
    HeliumFile();
    CarbonFile();
    OxygenFile();   
    NeonFile();
    MagnesiumFile();   
    SiliconFile();
    NickelFile();
    NuclearEnergyFile();

    Luminosity();

    CentralInfo();
	
    EnergyFile();
}


void OutDynamics::AllVariables(const double &t)
{
    OutRadius << t << "  ";
    OutVelocity << t << "  ";
    OutMass << t << "  ";
    OutDensity << t << "  ";
    OutTemperature << t << "  ";
    OutPressure << t << "  ";
    OutInternalEnergy << t << "  ";
    OutSoundSpeed << t << "  ";

    OutHelium << t << tab;
    OutCarbon << t << tab;
    OutOxygen << t << tab;
    OutNeon << t << tab;
    OutMagnesium << t << tab;
    OutSilicon << t << tab;
    OutNickel << t << tab;

    OutNuclearEnergy << t << tab;

    for (std::vector<int>::size_type j = 0; j < star_ptr->Nshells; j++){
        OutRadius << star_ptr->shell[j].Rout << "  ";
        OutVelocity << star_ptr->shell[j].Rdot_out << "  ";
        OutMass << star_ptr->shell[j].M << "  ";
        OutDensity << star_ptr->shell[j].Rho << "  ";
        OutTemperature << star_ptr->shell[j].T << "  ";
        OutPressure << star_ptr->shell[j].P << "  ";
        OutInternalEnergy << star_ptr->shell[j].u << "  ";
        OutSoundSpeed << star_ptr->shell[j].cs << "  ";

        OutHelium << star_ptr->shell[j].X[0] << "  ";
        OutCarbon << star_ptr->shell[j].X[1] << "  ";
        OutOxygen << star_ptr->shell[j].X[2] << "  ";
        OutNeon << star_ptr->shell[j].X[3] << "  ";
        OutMagnesium << star_ptr->shell[j].X[4] << "  ";
        OutSilicon << star_ptr->shell[j].X[5] << "  ";
        OutNickel << star_ptr->shell[j].X[6] << "  ";
        OutNuclearEnergy << star_ptr->shell[j].Eps << "  ";

    }

    OutRadius << std::endl;
    OutVelocity << std::endl;
    OutMass << std::endl;
    OutDensity << std::endl;
    OutTemperature << std::endl;
    OutPressure << std::endl;
    OutInternalEnergy << std::endl;
    OutSoundSpeed << std::endl;

    OutHelium << std::endl;
    OutCarbon << std::endl;
    OutOxygen << std::endl;
    OutNeon << std::endl;
    OutMagnesium << std::endl;
    OutSilicon << std::endl;
    OutNickel << std::endl;
    OutNuclearEnergy << std::endl;

    OutEnergy << t << tab << star_ptr->Kinetic_Energy() << tab 
                << star_ptr->Gravitational_Energy() << tab 
                << star_ptr->Internal_Energy() << tab 
                << star_ptr->Nuclear_Energy() << tab 
                << (star_ptr->Gravitational_Energy() 
                    + star_ptr->Internal_Energy() 
                    + star_ptr->Kinetic_Energy()) << std::endl;

    OutCentralInfo << t << tab << log10(star_ptr->shell[0].Rho) << tab 
                << log10(star_ptr->shell[0].T) << std::endl;

    OutLuminosity << t << tab << star_ptr->Luminosity << std::endl;
}



void OutDynamics::RadiusFile()
{
    std::string FileName = star_ptr->OutputFilePath + "radius.dat";
    OutRadius.open (FileName.c_str());
    OutRadius << "# This file contains the informations about radii of system" << std::endl;
    OutRadius << "# Author: Elvis Soares" << std::endl;
    OutRadius << "# time(s)   Radius[i](cm)" << std::endl;
}

void OutDynamics::VelocityFile()
{
    std::string FileName = star_ptr->OutputFilePath + "velocity.dat";
    OutVelocity.open (FileName.c_str());
    OutVelocity << "# This file contains the radial velocity of each shell as a function of time" << std::endl;
    OutVelocity << "# Author: Elvis Soares" << std::endl;
    OutVelocity << "# time(s)   Velocity[i](cm/s)" << std::endl;
}

void OutDynamics::MassFile()
{
    std::string FileName = star_ptr->OutputFilePath + "mass.dat";
    OutMass.open (FileName.c_str());
    OutMass << "# This file contains the mass of each shell as a function of time" << std::endl;
    OutMass << "# Author: Elvis Soares" << std::endl;
    OutMass << "# time(s)   Mass[i](g)" << std::endl;
}

void OutDynamics::DensityFile()
{
    std::string FileName = star_ptr->OutputFilePath + "density.dat";
    OutDensity.open (FileName.c_str());
    OutDensity << "# This file contains the density of each shell as a function of time" << std::endl;
    OutDensity << "# Author: Elvis Soares" << std::endl;
    OutDensity << "# time(s)   Density[i](cm)" << std::endl;
}

void OutDynamics::TemperatureFile()
{
    std::string FileName = star_ptr->OutputFilePath + "temperature.dat";
    OutTemperature.open (FileName.c_str());
    OutTemperature << "# This file contains the temperature of each shell as a function of time" << std::endl;
    OutTemperature << "# Author: Elvis Soares" << std::endl;
    OutTemperature << "# time(s)   Temperature[i](K)" << std::endl;
}

void OutDynamics::PressureFile()
{
    std::string FileName = star_ptr->OutputFilePath + "pressure.dat";
    OutPressure.open (FileName.c_str());
    OutPressure << "# This file contains the pressure of each shell as a function of time" << std::endl;
    OutPressure << "# Author: Elvis Soares" << std::endl;
    OutPressure << "# time(s)   Pressure [i](erg/cm³)" << std::endl;
}

void OutDynamics::InternalEnergyFile()
{
    std::string FileName = star_ptr->OutputFilePath + "internalenergy.dat";
    OutInternalEnergy.open (FileName.c_str());
    OutInternalEnergy << "# This file contains the internal energy of each shell as a function of time" << std::endl;
    OutInternalEnergy << "# Author: Elvis Soares" << std::endl;
    OutInternalEnergy << "# time(s)   u [i](erg/g)" << std::endl;
}

void OutDynamics::SoundSpeedFile()
{
    std::string FileName = star_ptr->OutputFilePath + "soundspeed.dat";
    OutSoundSpeed.open (FileName.c_str());
    OutSoundSpeed << "# This file contains the pressure of each shell as a function of time" << std::endl;
    OutSoundSpeed << "# Author: Elvis Soares" << std::endl;
    OutSoundSpeed << "# time(s)   Sound Speed [i](cm/s)" << std::endl;
}

void OutDynamics::HeliumFile()
{
    std::string FileName = star_ptr->OutputFilePath + "helium.dat";
    OutHelium.open (FileName.c_str());
    OutHelium << "# This file contains the helium abundance of each shell as a function of time" << std::endl;
    OutHelium << "# Author: Elvis Soares" << std::endl;
    OutHelium << "# time(s)   Abundance[i]( )" << std::endl;
}

void OutDynamics::CarbonFile()
{
    std::string FileName = star_ptr->OutputFilePath + "carbon.dat";
    OutCarbon.open (FileName.c_str());
    OutCarbon << "# This file contains the carbon abundance of each shell as a function of time" << std::endl;
    OutCarbon << "# Author: Elvis Soares" << std::endl;
    OutCarbon << "# time(s)   Abundance[i]( )" << std::endl;
}

void OutDynamics::OxygenFile()
{
    std::string FileName = star_ptr->OutputFilePath + "oxygen.dat";
    OutOxygen.open (FileName.c_str());
    OutOxygen << "# This file contains the oxygen abundance of each shell as a function of time" << std::endl;
    OutOxygen << "# Author: Elvis Soares" << std::endl;
    OutOxygen << "# time(s)   Abundance[i]( )" << std::endl;
}

void OutDynamics::NeonFile()
{
    std::string FileName = star_ptr->OutputFilePath + "neon.dat";
    OutNeon.open (FileName.c_str());
    OutNeon << "# This file contains the neon abundance of each shell as a function of time" << std::endl;
    OutNeon << "# Author: Elvis Soares" << std::endl;
    OutNeon << "# time(s)   Abundance[i]( )" << std::endl;
}

void OutDynamics::MagnesiumFile()
{
    std::string FileName = star_ptr->OutputFilePath + "magnesium.dat";
    OutMagnesium.open (FileName.c_str());
    OutMagnesium << "# This file contains the magnesium abundance of each shell as a function of time" << std::endl;
    OutMagnesium << "# Author: Elvis Soares" << std::endl;
    OutMagnesium << "# time(s)   Abundance[i]( )" << std::endl;
}

void OutDynamics::NickelFile()
{
    std::string FileName = star_ptr->OutputFilePath + "nickel.dat";
    OutNickel.open (FileName.c_str());
    OutNickel << "# This file contains the nickel abundance of each shell as a function of time" << std::endl;
    OutNickel << "# Author: Elvis Soares" << std::endl;
    OutNickel << "# time(s)   Abundance[i]( )" << std::endl;
}

void OutDynamics::SiliconFile()
{
    std::string FileName = star_ptr->OutputFilePath + "silicon.dat";
    OutSilicon.open (FileName.c_str());
    OutSilicon << "# This file contains the silicon abundance of each shell as a function of time" << std::endl;
    OutSilicon << "# Author: Elvis Soares" << std::endl;
    OutSilicon << "# time(s)   Silicon Abundance[i]( )" << std::endl;
}

void OutDynamics::NuclearEnergyFile()
{
    std::string FileName = star_ptr->OutputFilePath + "nuclearenergy.dat";
    OutNuclearEnergy.open (FileName.c_str());
    OutNuclearEnergy << "# This file contains the nuclear energy rate of each shell as a function of time" << std::endl;
    OutNuclearEnergy << "# Author: Elvis Soares" << std::endl;
    OutNuclearEnergy << "# time(s)   Eps(erg/g s)" << std::endl;
}

void OutDynamics::EnergyFile()
{
    std::string FileName = star_ptr->OutputFilePath + "energy.dat";
    OutEnergy.open (FileName.c_str());
    OutEnergy << "# This file contains the informations about energy of system" << std::endl;
    OutEnergy << "# Author: Elvis Soares" << std::endl;
    OutEnergy << "# time(s)   Kinectic (erg)    Gravitational(erg)     Internal(erg)    Nuclear_Energy(erg)    Total Energy(erg)" << std::endl;
}

void OutDynamics::Luminosity()
{
    std::string FileName = star_ptr->OutputFilePath + "luminosity.dat";
    OutLuminosity.open (FileName.c_str());
    OutLuminosity << "# This file contains the informations about the superficial Luminosity of Supernova" << std::endl;
    OutLuminosity << "# Author: Elvis Soares" << std::endl;
    OutLuminosity << "# time(s)   L (erg/s)" << std::endl;
}

void OutDynamics::CentralInfo()
{
    std::string FileName = star_ptr->OutputFilePath + "centralinfo.dat";
    OutCentralInfo.open (FileName.c_str());
    OutCentralInfo << "# This file contains the informations about the centre of Supernova" << std::endl;
    OutCentralInfo << "# Author: Elvis Soares" << std::endl;
    OutCentralInfo << "# time(s)   log Rho (g/cm³)    log T (K)" << std::endl;
}
