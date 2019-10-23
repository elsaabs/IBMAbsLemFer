#ifndef POPULATION_H
#define POPULATION_H

#include <iostream>
#include <fstream>      /* pour ecrire dans un fichier */
#include <cmath>        /* outils mathematiques */
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
        Population(int const& squareLength, double const& sizeInitialEnzyme, double const& sizeInitialC,
                 double const& sizeInitialD,
                 double const& indDeathZRate, double const& nbrEquivZD, double const& energCostZD,
                 double const& indDeathCRate, double const& creationCRate, double const& indBindZCRate,
                 double const& nbrEquivCD, double const& indDeathDRate, double const& creationDRate,
                 double const& erosion, double const& paramH, double const& diffusionD,
                 int const& nbrParCase, Microbe const& microbeInitial);


        void change(double const& sizeInitialEnzyme, double const& sizeInitialC, double const& sizeInitialD,
                 double const& sizeInitialEnzyme2, double const& sizeInitialC2, double const& sizeInitialD2,
                 int const& nbrParCase, int const& nbrParCase2, int const& numeroCaseInitiale, int const& nbrCasesMutant,
                 Microbe const& microbe1Initial, Microbe const& microbe2Initial);

        ~Population();

        //Z,C,D
        void oneStepEulerOneCase(double const& sizeStep, int const& index);
        void eulerSchemeOneCase(double const& sizeStep, double const& intervalTime, int const& index);
        void eulerScheme(double const& sizeStep, double const& intervalTime, std::vector<double> &birthTab);
        void diffusionCarbon(double const& sizeTimeStep);
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
        int ifMovingOneCase(int const& indexInitial, double const& phi, double const& probaReplace);
        void addMicrobe(double const& initialGrowth, int const& indexSubPopM, int const& position,
                            double const& probaMoving, double const& probaReplace);
        void deleteMicrobe(double &carbonCReleased, double &carbonDReleased, double const& indexMicrobe);
        void deleteMicrobe(double &carbonCReleased, double &carbonDReleased, int const& indexSubPop,
                               int const& indexInPop);
        void oneStepDiffusion(double &timeStep, double const& timeStepEuler, double const& probaMoving,
                              double const& diffusionTimeStep, double const& probaReplace);
        double getphi(int const& indexSubPop) const;
        void calculusStepMeanVar(double const& timeStep, double &meanZ, double &varZ,
                                 double &meanC, double &varC, double &meanD,
                                 double &varD, double &meanM, double &varM,
                                 double &meanaccessDOC, double &varaccessDOC) const;
        void calculusMeanVarBiomass(double &meanZ, double &varZ,
                                 double &meanC, double &varC, double &meanD,
                                 double &varD, double &meanM, double &varM,
                                 double &meanaccessDOC, double &varaccessDOC,
                                 double const& normalParamK, double const& weightEnzyme,
                                 double const& weightCarbonC, double const& weightCarbonD) const;


        void writeFile(std::string nameFile, double const& timeIndex, double const& normalParamK,
                       double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD);
        void writeFile2types(std::string nameFile, double const& timeIndex, double const& normalParamK,
                           double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD,
                           double const& phi1, double const& phi2);
        void writeFile2typesTotalAccess(std::string nameFile, double const& timeIndex, double const& normalParamK,
                           double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD,
                           double const& phi1, double const& phi2);
        void writeTotalAccess(std::string nameFile, double const& timeIndex, double const& normalParamK,
                           double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD,
                                  double const& write1, double const& write2, double const& write3, double const& write4, 
                                  double const& write5, double const& write6, double const& write7, double const& write8, 
                                  double const& write9, double const& write10);






    private:
        //parametres pour C,Z,D
        int m_length; //square side size
        std::vector<double> m_sizeZ; //total size of enzyme species in each box
        std::vector<double> m_sizeC; //total size of carbonC species in each box
        std::vector<double> m_sizeD; //total size of carbonD species in each box
        std::vector<int> m_sizeM;
        std::vector<double> m_linearGrowthtempInd; //individual linear growth with Euler scheme
        std::vector<double> m_types; //gives the types in each case
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
        double m_diffusionD; //gives the diffusion coefficient of entity D
        std::vector<double> m_productionZ; //gives the production rate by the pop of microbes, careful: initialized to -1!!
        std::vector<double> m_maxlinearGrowth; //gives the individual growth rate in each box, careful: initialized to -1!!
        double m_gammaZ; //attention, ici on suppose que tous les ind ont le meme gammaZ (programme moins lourd...)
        double m_maxCostD; //attention, ici on suppose que tous les ind ont le meme maxCostD (programme moins lourd...)
        double m_maxGrowthrateforanyphi; //attention, ici on suppose que tous les ind ont le meme maxCostD (programme moins lourd...)
        //std::vector<double> m_consumptionDprodZ; //gives the consumption rate of D by the pop of microbes to create Z, careful: initialized to -1!!
        //std::vector<double> m_consumptionDprodM;

        //parametres pour M
        std::vector<std::vector<double> > m_individuals; //represents the pop of microbes in the model,
                                                 //each subpopulation has identical parameters
                                                 // each subPop is formed by a vector with the quantity of stocked D by the microbe
        std::vector<std::vector<int> > m_positions; //and positions
        std::vector<Microbe> m_microbeTypes; //microbe type for each column of the previous table.
        double m_sizeMTotal; //retient le nombre total d'individu M (evite de le recalculer a chaque pas de temps...)
        double m_maxdeathrate; //retient le taux INDIVIDUEL max
        //double m_maxabsorptionrate; //retient le taux INDIVIDUEL max
        double m_timeStepDiffusion; //pas de temps pour la diffusion

};



#endif // POPULATION_H

