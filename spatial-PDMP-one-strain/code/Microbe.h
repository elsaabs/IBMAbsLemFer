#ifndef MICROBE_H
#define MICROBE_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <assert.h>



class Microbe
{
    public:
        Microbe(double const& nbrEquivCarbonD, double const& nbrLooseCarbonD,
               double const& phi, double const& maxIndBirthRate, double const& paramH,
               double const& indDeathRate, double const& gammaM,
               double const& gammaZ, double const& paramP);
        Microbe(Microbe const& microbeCopy, double const& newPhi);
        ~Microbe();

        double deathRate() const;
        double maxdeathRate() const;
        double maxproduceZRate() const;
        double maxlineargrowthM() const;
        void carbonReleasedwhenDeath(double &carbonCReleased, double &carbonDReleased,
                                     double const& equivCD, double const& actualGrowth) const;
        double getphi() const;
        double getgammaZ() const;
        double getmaxCost() const;
        double getmaxGrowthforanyphi() const;


    private:
        double m_nbrEquivCarbonD; //correspond to "eta" -> structural cost (eta=alpha/K)
        double m_nbrLooseCarbonD; //correspond to "eta'" -> energetic cost (eta'=alpha'/K)
        double m_phi; //repartition reproduction/production of Z
        double m_maxIndBirthRate; //correspond to parameter "a"
        double m_paramH; //correspond to parameter "h"
        double m_indDeathRate; //correspond to parameter "d_M"
        double m_gammaM; // = 1/(alpha+alpha') -> probability to get a microbe during a respiration
        double m_gammaZ; // = 1/(beta+beta') -> probability to get an enzyme during a respiration
        double m_paramP; //parameter p -> proportion of carbonC w.r.t. carbonD
        int m_nbrMaxEventAbs; //give the number of events of absorption before getting a new individual


};


#endif // MICROBE_H
