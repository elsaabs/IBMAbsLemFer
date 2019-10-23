#include <iostream>
#include "./code/Population.h"
#include <string>
#include <sstream>
#include <math.h>

using namespace std;

int main(int argc, char* argv[])
{


    //Argument
    cout << "argc = " << argc << endl;
    if (argc != 19)
    {
        cout << "18 arguments dans l'ordre: 1-ID cluster, 2-ID job, 3-timeMax, 4- initalTimeforMeans, 5-phi,";
        cout << "6-diffusionD, 7-VmU, 8-KmU, 9-dM, 10-dZ, 11-gammaM, 12-gammaZ,";
        cout << " 13-paramP, 14-VmD, 15-InputC, 16-InputD, 17-LeachingC, 18-LeachingD" << endl;
        return -1;
    }
    
    int const squareLength(10);
    double const diffusionTimeStep(1e-3);
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //parameters to specify
    //initialisation of the stream
    srand(atof(argv[1])+atof(argv[2]));
    
    
    //parametre de normalisation K
    //L'idee de ce programme est de pouvoir regarder des K grands typiquement K~alpha=wM/wD, donc ne pas regarder des micro bouts d'espace...
    //Ce qu'on ne peut pas faire avec le programme 4esp car l'algo est trop lourd!
    double const normalParamK(10); // -> nbr d'individus moyen par case
    
    //parameters of the entities:
    double const erosion(0);
    //parameters of cost:
    double const weightMicrobe(1e-9);  //wM -> biomasse d'1 microbe (i.e. poids d'1 microbe, dans la mesure de poids que tu souhaites obtenir a la fin: ex: mg) = 10^-9 mg
    double const weightCarbonC(1e-16);  //wC -> biomasse d'1 carbonC (i.e. poids d'1 C) = 10^-16 mg
    double const weightEnzyme(1e-16);  //wZ -> biomasse d'1 enzyme (i.e. poids d'1 Z) = 10^-16 mg
    double const weightCarbonD(1e-19);  //wD -> biomasse d'1 carbonD (i.e. poids d'1 D)
    //on doit avoir alpha=wM/wD > beta=wC/wD > gamma=wZ/wD
    
    //Microbe
    double const VmU(atof(argv[7]));     //(0.42);
    double const KmU(atof(argv[8]));
    double const indDeathRateM(atof(argv[9]));     //(2e-4); //dM
    double const gammaM(atof(argv[11]));   //(0.3);
    double const gammaZ(atof(argv[12]));
    double const paramP(atof(argv[13]));  //p
    double const phi(atof(argv[5]));  //phi initial...
    
    //Enzyme
    double const indDeathRateZ(atof(argv[10]));  //dZ
    //double biomasseInitialZ(1.4e-10);  // -> en biomasse totale par case
    
    
    //CarbonC
    double const IC(atof(argv[15]));
    double const indDeathCRate(atof(argv[17]));  //l_C
    double const VmD(atof(argv[14]));
    //double biomasseInitialC(1.2e-7);  // -> en biomasse totale par case
    
    
    //CarbonD
    double const ID(atof(argv[16]));
    double const indDeathDRate(atof(argv[18]));  //l_D
    double const diffusionD(atof(argv[6])); //en cm^2.h-1
    double const diffusionDnormal(diffusionD*1e6/pow(normalParamK,0.666)); //ne pas toucher
    
    //double biomasseInitialD(3.0e-12);  //  -> en biomasse totale par case
    
    //mouvement
    double const probaMoving(0.3); //probability to move at the birth time
    double const probaReplace(0.01);
    
    //Time and file parameters to be adjusted
    double const nbrWritting(1000);  //Numbers of intervals time to write on the file
    double const timeMax(atof(argv[3]));  //Maximal time to observe
    double const timeStepEuler(0.001);  //pas de resolution des EDO


    ostringstream oss;
    oss << "./cases/scoresCases_" << argv[1] << "_" << argv[2] << ".txt";
    string nameFile(oss.str());   //output namefile
    ofstream flux(nameFile.c_str());

    ostringstream oss1;
    oss1 << "./totaux/scoresTotaux_" << argv[1] << "_" << argv[2] << ".txt";
    string nameFileTotal(oss1.str());   //output namefile
    ofstream fluxTotal(nameFileTotal.c_str());    //on met cette commande pour effacer les anciennes
                                                  //donnees dans les fichiers, s'il y en a!

    //end: parameters to specify
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    //others parameters
    double const structuralCostCarbonDMicrobeNormal((weightMicrobe/ weightCarbonD) / normalParamK );
    double const energeticCostCarbonDMicrobeNormal((1-gammaM)/gammaM * structuralCostCarbonDMicrobeNormal);
    double const structuralCostCarbonDCarbonC(weightCarbonC/weightCarbonD);
    double const structuralCostCarbonDEnzyme(weightEnzyme/weightCarbonD);
    double const energeticCostCarbonDEnzyme((1-gammaZ)/gammaZ * structuralCostCarbonDEnzyme);
    double const bargammaZ(gammaZ / structuralCostCarbonDEnzyme);
    double const creationCRate(IC / weightCarbonC);
    double const creationDRate(ID / weightCarbonD);
    double const indBindCZRate(VmD * weightEnzyme);
    double const paramH(KmU / weightCarbonD);
    //les 4 valeurs suivantes seront initialisees par le membre eqAttendus.
    int nbrInitialMParCase(0);  // Initial number of Microbes -> en nbr d'individus total par case
    double sizeInitialZ(0); //(biomasseInitialZ/(normalParamK*weightEnzyme));
    double sizeInitialC(0); //(biomasseInitialC/(normalParamK*weightCarbonC));
    double sizeInitialD(0); //(biomasseInitialD/(normalParamK*weightCarbonD));




    //Microb(structuralCostM,energeticCostM,phi,maxIndBirthRate(VmU),paramH(KmU),indDeathRate(dm),gammaM,gammaZ,p)
    Microbe microbe1Init(structuralCostCarbonDMicrobeNormal,energeticCostCarbonDMicrobeNormal,phi,VmU,paramH,indDeathRateM,gammaM,
                        bargammaZ,paramP);

    eqAttendus(phi, gammaM, VmU, indDeathRateM, gammaZ, indDeathRateZ, KmU, VmD, indDeathCRate, indDeathDRate,
                 erosion, IC, ID, paramP, weightMicrobe, weightEnzyme, weightCarbonC, weightCarbonD, normalParamK,
                 nbrInitialMParCase, sizeInitialZ, sizeInitialC, sizeInitialD);

    //construction 2 -> un type de bacteries identiques partout sur la grille
    Population popTotal(squareLength,sizeInitialZ,sizeInitialC,sizeInitialD,indDeathRateZ,structuralCostCarbonDEnzyme,
                        energeticCostCarbonDEnzyme,indDeathCRate,creationCRate,indBindCZRate,structuralCostCarbonDCarbonC,
                        indDeathDRate,creationDRate,erosion,paramH,diffusionDnormal,nbrInitialMParCase,microbe1Init);



    //other time and file parameters
    double timeIndex(0);
    double timeStep(0);
    bool existence(true);
    double const initialTimeForMeans(atof(argv[4]));
    if (initialTimeForMeans >= timeMax) {cout << "probleme dans l'initialisation de initialTimeForMeans"; abort();}
    double meanZ(0);
    double varZ(0);
    double meanC(0);
    double varC(0);
    double meanD(0);
    double varD(0);
    double meanM(0);
    double varM(0);
    double meanaccessDOC(0);
    double varaccessDOC(0);
    
    double meanZ2(0);
    double varZ2(0);
    double meanC2(0);
    double varC2(0);
    double meanD2(0);
    double varD2(0);
    double meanM2(0);
    double varM2(0);
    double meanaccessDOC2(0);
    double varaccessDOC2(0);

    popTotal.calculusMeanVarBiomass(meanZ2,varZ2,meanC2,varC2,meanD2,varD2,meanM2,varM2,meanaccessDOC2,varaccessDOC2,
                                        normalParamK,weightEnzyme,weightCarbonC,weightCarbonD);

    popTotal.writeFile(nameFile,timeIndex,normalParamK,weightEnzyme,weightCarbonC,weightCarbonD);
    popTotal.writeTotalAccess(nameFileTotal,timeIndex,normalParamK,weightEnzyme,weightCarbonC,weightCarbonD,
                              meanZ2,varZ2,meanC2,varC2,meanD2,varD2,meanM2,varM2,meanaccessDOC2,varaccessDOC2);


    for (int k(1); k <= nbrWritting; k++)
    {
        while (timeIndex < (timeMax / nbrWritting)*k && existence)
        {
            popTotal.oneStepDiffusion(timeStep,timeStepEuler,probaMoving,diffusionTimeStep,probaReplace);
            timeIndex += timeStep;

            //calcul moyenne et variance pour les sorties de communs
            if (timeIndex > initialTimeForMeans)
            {
                popTotal.calculusStepMeanVar(timeStep,meanZ,varZ,meanC,varC,meanD,varD,meanM,varM,meanaccessDOC,varaccessDOC);
            }


            if (popTotal.sizePopM() == 0) //plus de bacteries
            {
                cout << "stop: no bacteria after time " << timeIndex << endl;
                existence = false;
                k = nbrWritting;
            }
            if (popTotal.sizePopM() > normalParamK * 10000) //population trop grande: probleme dans la programmation
            {
                cout << "stop: too much microbes! After time " << timeIndex << endl;
                existence = false;
                k = nbrWritting;
            }

        }

        popTotal.calculusMeanVarBiomass(meanZ2,varZ2,meanC2,varC2,meanD2,varD2,meanM2,varM2,meanaccessDOC2,varaccessDOC2,
                                        normalParamK,weightEnzyme,weightCarbonC,weightCarbonD);

        popTotal.writeFile(nameFile,timeIndex,normalParamK,weightEnzyme,weightCarbonC,weightCarbonD);
        popTotal.writeTotalAccess(nameFileTotal,timeIndex,normalParamK,weightEnzyme,weightCarbonC,weightCarbonD,
                                  meanZ2,varZ2,meanC2,varC2,meanD2,varD2,meanM2,varM2,meanaccessDOC2,varaccessDOC2);
    }

    double popTotalFinal((double)popTotal.sizePopM());
    meanZ = meanZ*normalParamK*weightEnzyme/(timeMax-initialTimeForMeans);
    varZ = varZ*normalParamK*normalParamK*weightEnzyme*weightEnzyme/(timeMax-initialTimeForMeans)-meanZ*meanZ;
    meanC = meanC*normalParamK*weightCarbonC/(timeMax-initialTimeForMeans);
    varC = varC*normalParamK*normalParamK*weightCarbonC*weightCarbonC/(timeMax-initialTimeForMeans)-meanC*meanC;
    meanD = meanD*normalParamK*weightCarbonD/(timeMax-initialTimeForMeans);
    varD = varD*normalParamK*normalParamK*weightCarbonD*weightCarbonD/(timeMax-initialTimeForMeans)-meanD*meanD;
    meanM = meanM/(timeMax-initialTimeForMeans);
    varM = varM/(timeMax-initialTimeForMeans)-meanM*meanM;
    meanaccessDOC = meanaccessDOC*normalParamK*weightCarbonD/(timeMax-initialTimeForMeans);
    varaccessDOC = varaccessDOC*normalParamK*normalParamK*weightCarbonD*weightCarbonD/(timeMax-initialTimeForMeans)
                   -meanaccessDOC*meanaccessDOC;

    ostringstream oss2;
    oss2 << "./communs/scoresCommuns" << argv[1] << "_" << argv[2] << ".txt";
    string nameFileCommun(oss2.str());
    ofstream fluxcommun(nameFileCommun.c_str()); //efface les données présentes (si anciennes simus)
    fluxcommun << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5];
    fluxcommun << " " << argv[6] << " " << argv[7] << " " << argv[8] << " " << argv[9] << " " << argv[10];
    fluxcommun << " " << argv[11] << " " << argv[12] << " " << argv[13] << " " << argv[14] << " " << argv[15];
    fluxcommun << " " << argv[16] << " " << argv[17] << " " << argv[18];
    fluxcommun << " " << timeIndex << " " << popTotalFinal << " "  << meanZ  << " " << VmD*meanZ << " " << varZ;
    fluxcommun << " " << meanC << " " << varC << " " << meanD << " " << varD << " ";
    fluxcommun << meanM << " " << varM << " " << meanaccessDOC << " " << varaccessDOC << endl;

    return 0;
}
