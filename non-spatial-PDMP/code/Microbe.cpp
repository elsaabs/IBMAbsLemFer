#include "Microbe.h"


using namespace std;

Microbe::Microbe(double const& nbrEquivCarbonD, double const& nbrLooseCarbonD,
               double const& phi, double const& maxIndBirthRate, double const& paramH,
               double const& indDeathRate, double const& gammaM,
               double const& gammaZ, double const& paramP):
               m_nbrEquivCarbonD(nbrEquivCarbonD), m_nbrLooseCarbonD(nbrLooseCarbonD),
               m_maxIndBirthRate(maxIndBirthRate), m_paramH(paramH),
               m_indDeathRate(indDeathRate), m_paramP(paramP)
{
    if ((phi <= 1)&&(phi >= 0)&&(gammaM <= 1)&&(gammaM >= 0)&&(gammaZ <= 1)&&(gammaZ >= 0)&&(paramP <= 1)&&(paramP >= 0))
    {
        m_phi = phi;
        m_gammaM = gammaM;
        m_gammaZ = gammaZ;
        m_paramP = paramP;
    }

    else
    {
        cout << "problem with the definition of phi, gammaZ, gammaM or paramP";
        assert(false);
    }
}






Microbe::Microbe(Microbe const& microbeCopy, double const& newPhi)
{
    m_nbrEquivCarbonD = microbeCopy.m_nbrEquivCarbonD;
    m_nbrLooseCarbonD = microbeCopy.m_nbrLooseCarbonD;
    m_phi = newPhi;
    m_maxIndBirthRate = microbeCopy.m_maxIndBirthRate;
    m_paramH = microbeCopy.m_paramH;
    m_indDeathRate = microbeCopy.m_indDeathRate;
    m_gammaM = microbeCopy.m_gammaM;
    m_gammaZ = microbeCopy.m_gammaZ;
    m_paramP = microbeCopy.m_paramP;
}



Microbe::~Microbe()
{
    //dtor
}





double Microbe::deathRate() const
{
    return m_indDeathRate;
}


double Microbe::maxdeathRate() const
{
    return m_indDeathRate;
}



double Microbe::maxproduceZRate() const
{
    return m_phi * m_gammaZ * m_nbrEquivCarbonD * m_maxIndBirthRate;
}


double Microbe::maxlineargrowthM() const
{
    return (1-m_phi) * m_gammaM * m_maxIndBirthRate;
}





void Microbe::carbonReleasedwhenDeath(double &carbonCReleased, double &carbonDReleased, double const& equivCD, double const& linearGrowth) const
{
    carbonDReleased = (1-m_paramP) * m_nbrEquivCarbonD * (1 + linearGrowth);
    carbonCReleased = (m_paramP) * m_nbrEquivCarbonD / equivCD * (1 + linearGrowth);
}



double Microbe::getphi() const
{
    return m_phi;
}


double Microbe::getgammaZ() const
{
    return m_gammaZ;
}


double Microbe::getmaxCost() const
{
    return m_nbrEquivCarbonD+m_nbrLooseCarbonD;
}


double Microbe::getmaxGrowthforanyphi() const
{
    return m_gammaM * m_maxIndBirthRate;
}



