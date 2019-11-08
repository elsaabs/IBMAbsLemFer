#ifndef POPULATION_H
#define POPULATION_H

#include <iostream>
#include <fstream>      /* to write in a file */
#include <cmath>        /* mathematical tools */
#include <algorithm>
#include <cstdlib>      /* abort */
#include <stdlib.h>     /* srand, rand */
#include <vector>
#include <assert.h>
#include "Microbe.h"


int varUnif(int const& maxInt);
double varUniform();
double varExpo(double const& parameter);
void eqAttendus(double const& phi, double const& gammaM, double const& VmU, double const& dM,
                  double const& gammaZ, double const& dZ, double const& KmU, double const& VmD, double const& lC, double const& lD,
                  double const& erosion, double const& IC, double const& ID, double const& paramP, double const& weightMicrobe,
                  double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD, double const& normalParamK,
                  int &MeqInd, double &ZeqDens, double &CeqDens, double &DeqDens);



class Population
{
    public:
        Population(double const& sizeInitialEnzyme, double const& sizeInitialC,
                 double const& sizeInitialD, double const& indDeathZRate, double const& nbrEquivZD, double const& energCostZD,
                 double const& indDeathCRate, double const& creationCRate, double const& indBindZCRate,
                 double const& nbrEquivCD, double const& indDeathDRate, double const& creationDRate,
                 double const& erosion, double const& paramH, int const& nbrInitialMicrobes, Microbe const& microbeInitial);
        ~Population();

        //Z,C,D
        void oneStepEuler(double const& sizeStep);
        void eulerScheme(double const& sizeStep, double const& intervalTime, std::vector<double> &birthTab);
        double totalZ() const;
        double totalC() const;
        double totalD() const;
        //M
        int nbrTypes() const;
        double sizePopM() const;
        double sizeSubPopM(int const& indexType) const;
        void functionIndexM(double const& indexMicrobe, int &indexSubPopM, int &indexInPopM) const;
        double maxDeathRate();
        double deathRate(double const& indexMicrobe) const;
        void addMaxRate(int const& indexSubPopM);
        void suppMaxRate(double const& death);
        void updateTotals(int const& indexSubPopM, int const& indexInPopM, bool const& isplus);
        void addMicrobe(double const& initialGrowth, int const& indexSubPopM, int const& indexInPopM,
                            double const& probaMutation, double const& paramMutation);
        void deleteMicrobe(double &carbonCReleased, double &carbonDReleased, double const& indexMicrobe);
        void oneStepAlgo(double &timeStep, double const& timeStepEuler, double const& probaMutation,
                                    double const& paramMutation, double const& annexTimeStep);
        double getphi(int const& indexSubPop) const;
        void writeFile(std::string nameFile, double const& timeIndex, double const& normalParamK,
                       double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD);
        void writeFile2types(std::string nameFile, double const& timeIndex, double const& normalParamK,
                           double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD,
                           double const& phi1, double const& phi2);
     






    private:
        //parametres for C,Z,D
        double m_sizeZ; //total size of enzyme species
        double m_sizeC; //total size of carbonC species
        double m_sizeD; //total size of carbonD species
        double m_linearGrowthtempInd; //individual linear growth with Euler scheme
        double m_indDeathZRate; //correspond to d_Z
        double m_nbrEquivZD; //correspond to "gamma" -> structural cost
        double m_energCostZD; //correspond to "gamma'" -> energetic cost.
        double m_indDeathCRate; //correspond to l_C
        double m_creationCRate; //correspond to I_C
        double m_indBindZCRate; //correspond to the parameter before zc
        double m_nbrEquivCD; //correspond to "beta" -> structural cost
        double m_indDeathDRate; //correspond to l_D
        double m_creationDRate; //correspond to I_D
        double m_erosion; //correspond to epsilon
        double m_paramH; //gives the normalized parameter KmU.
        double m_productionZ; //gives the production rate by the pop of microbes, careful: initialized to -1!!
        double m_maxlinearGrowth; //gives the individual growth rate, careful: initialized to -1!!
        double m_gammaZ; //we assume that all individuals have the same gammaZ
        double m_maxCostD; //we assume that all individuals have the same maxCostD
        double m_maxGrowthrateforanyphi; //we assume that all individuals have the same maxCostD

        //parametres for M
        std::vector<std::vector<double> > m_individuals; //represents the pop of microbes in the model,
                                                 //each subpopulation has identical parameters
                                                 // each subPop is formed by a vector with the quantity of stocked D by the microbe
        std::vector<Microbe> m_microbeTypes; //microbe type for each column of the previous table.
        double m_sizeMTotal; //save the total number of individuals M (so that we don't have to recalculate it at each time step)
        double m_maxdeathrate; //save the maximum INDIVIDUAL rate
        double m_timeStepDiffusion; //time step of diffusion

};



#endif // POPULATION_H

