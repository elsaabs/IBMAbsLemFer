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
            cout << "18 arguments: 1-ID cluster, 2-ID job, 3-timeMax, 4-phi1 (resident),";
            cout << " 5-phi2 (mutant), 6-diffusionD, 7-VmU, 8-KmU, 9-dM, 10-dZ, 11-gammaM, 12-gammaZ,";
            cout << " 13-paramP, 14-VmD, 15-InputC, 16-InputD, 17-LeachingC, 18-LeachingD" << endl;
            return -1;
        }
    
    int const squareLength(10);
    double const diffusionTimeStep(1e-3);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //parameters to specify
    //initialisation of the stream
    srand(atof(argv[1])+atof(argv[2]));


    //parameter of normalisation K
    double const normalParamK(10); // -> average number of individuals per microsite

    //parameters of the entities:
    double const erosion(0);
    //parameters of cost:
    double const weightMicrobe(1e-9);  //wM -> biomass of 1 microbe (mg)
    double const weightCarbonC(1e-16); //wC -> biomass of 1 carbonC (mg) = typically cellulose
    double const weightEnzyme(1e-16);  //wZ -> biomass of 1 enzyme (mg) = typically cellulase
    double const weightCarbonD(1e-19);  //wD -> biomasse of 1 carbonD (mg) = typically glucose
    //we need to have alpha=wM/wD > beta=wC/wD > gamma=wZ/wD

    //Microbe
    double const VmU(atof(argv[7]));     //(0.42);
    double const KmU(atof(argv[8]));
    double const indDeathRateM(atof(argv[9]));     //(2e-4); //dM
    double const gammaM(atof(argv[11]));   //(0.3);
    double const gammaZ(atof(argv[12]));
    double const paramP(atof(argv[13]));  //p
    double const phi1(atof(argv[4]));  //value of resident phi
    double const phi2(atof(argv[5]));  //value of mutant phi
    int const nbrInitialCasesMutants(5);
    int const numeroCaseInitiale(((squareLength-1)/2)*squareLength+((squareLength-1)/2)); //number of initial microsite wheree we put mutants
    //int = natural number, double = real number

    //Enzyme
    double const indDeathRateZ(atof(argv[10]));  //dZ


    //CarbonC
    double const IC(atof(argv[15]));
    double const indDeathCRate(atof(argv[17]));  //l_C
    double const VmD(atof(argv[14]));


    //CarbonD
    double const ID(atof(argv[16]));
    double const indDeathDRate(atof(argv[18]));  //l_D
    double const diffusionD(atof(argv[6])); //en cm^2.h-1
    double const diffusionDnormal(diffusionD*1e6/pow(normalParamK,0.666)); //we won't change this value

    //movement
    double const probaMoving(0.3); //probability to move at the birth time
    double const probaReplace(0.01);

    //Time and file parameters to be adjusted
    double const nbrWritting(1000);  //Numbers of intervals time to write on the file
    double const timeMax(atof(argv[3]));  //Maximal time to observe
    double const timeStepEuler(0.001);  //step for solving differential equations

    //Parameter for initialisation
    int choix(1);

    ostringstream oss;
    oss << "./cases/scoresCases_" << argv[1] << "_" << argv[2] << ".txt";
    string nameFile(oss.str());   //output namefile
    ofstream flux(nameFile.c_str()); //to erase old data if there are some in the file

    ostringstream oss1;
    oss1 << "./totaux/scoresTotaux_" << argv[1] << "_" << argv[2] << ".txt";
    string nameFileTotal(oss1.str());   //output namefile
    ofstream fluxTotal(nameFileTotal.c_str());

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
    //the 4 following parameters are initiliased with "eqAttendus" (the deterministic steady state used as an estimation).
    int nbrInitialMParCase(0);  // Initial number of Microbes (total number per microsite)
    double sizeInitialZ(0); //(biomasseInitialZ/(normalParamK*weightEnzyme));
    double sizeInitialC(0); //(biomasseInitialC/(normalParamK*weightCarbonC));
    double sizeInitialD(0); //(biomasseInitialD/(normalParamK*weightCarbonD));


    bool stop(false);




    //Microb(structuralCostM,energeticCostM,phi,maxIndBirthRate(VmU),paramH(KmU),indDeathRate(dm),gammaM,gammaZ,p)
    Microbe microbe1Init(structuralCostCarbonDMicrobeNormal,energeticCostCarbonDMicrobeNormal,phi1,VmU,paramH,indDeathRateM,gammaM,
                        bargammaZ,paramP);

    Microbe microbe2Init(structuralCostCarbonDMicrobeNormal,energeticCostCarbonDMicrobeNormal,phi2,VmU,paramH,indDeathRateM,gammaM,
                        bargammaZ,paramP);

    if (choix == 1)
        {
            eqAttendus(phi1, gammaM, VmU, indDeathRateM, gammaZ, indDeathRateZ, KmU, VmD, indDeathCRate, indDeathDRate,
                 erosion, IC, ID, paramP, weightMicrobe, weightEnzyme, weightCarbonC, weightCarbonD, normalParamK,
                 nbrInitialMParCase, sizeInitialZ, sizeInitialC, sizeInitialD);
        }
    else if (choix == 2)
        {
            eqAttendus(phi2, gammaM, VmU, indDeathRateM, gammaZ, indDeathRateZ, KmU, VmD, indDeathCRate, indDeathDRate,
                 erosion, IC, ID, paramP, weightMicrobe, weightEnzyme, weightCarbonC, weightCarbonD, normalParamK,
                 nbrInitialMParCase, sizeInitialZ, sizeInitialC, sizeInitialD);
        }

    Population popTotal(squareLength,sizeInitialZ,sizeInitialC,sizeInitialD,indDeathRateZ,structuralCostCarbonDEnzyme,
                        energeticCostCarbonDEnzyme,indDeathCRate,creationCRate,indBindCZRate,structuralCostCarbonDCarbonC,
                        indDeathDRate,creationDRate,erosion,paramH,diffusionDnormal,nbrInitialMParCase,numeroCaseInitiale,
                        nbrInitialCasesMutants,microbe1Init,microbe2Init);

    if (choix == 3)
    {
        int nbrInitialMParCase2(0);
        double sizeInitialZ2(0);
        double sizeInitialC2(0);
        double sizeInitialD2(0);
        eqAttendus(phi1, gammaM, VmU, indDeathRateM, gammaZ, indDeathRateZ, KmU, VmD, indDeathCRate, indDeathDRate,
                 erosion, IC, ID, paramP, weightMicrobe, weightEnzyme, weightCarbonC, weightCarbonD, normalParamK,
                 nbrInitialMParCase, sizeInitialZ, sizeInitialC, sizeInitialD);
        eqAttendus(phi2, gammaM, VmU, indDeathRateM, gammaZ, indDeathRateZ, KmU, VmD, indDeathCRate, indDeathDRate,
                 erosion, IC, ID, paramP, weightMicrobe, weightEnzyme, weightCarbonC, weightCarbonD, normalParamK,
                 nbrInitialMParCase2, sizeInitialZ2, sizeInitialC2, sizeInitialD2);

        popTotal.change(sizeInitialZ,sizeInitialC,sizeInitialD,sizeInitialZ2,sizeInitialC2,
                            sizeInitialD2,nbrInitialMParCase,nbrInitialMParCase2,numeroCaseInitiale,
                        nbrInitialCasesMutants,microbe1Init,microbe2Init);
    }





    //other time and file parameters
    double timeIndex(0);
    double timeStep(0);
    double fraction(-1);
    bool coexistence(true);

    popTotal.writeFile2types(nameFile,timeIndex,normalParamK,weightEnzyme,weightCarbonC,weightCarbonD,phi1,phi2);
    popTotal.writeFile2typesTotalAccess(nameFileTotal,timeIndex,normalParamK,weightEnzyme,weightCarbonC,weightCarbonD,phi1,phi2);


    for (int k(1); k <= nbrWritting; k++)
    {
        while (timeIndex < (timeMax / nbrWritting)*k && coexistence)
        {
            popTotal.oneStepDiffusion(timeStep,timeStepEuler,probaMoving,diffusionTimeStep,probaReplace);
            timeIndex += timeStep;


            if (popTotal.nbrTypes() == 1) // if more than 1 type/strain of microbe
            {
                if (popTotal.getphi(0) == phi1)
                {
                    cout << "stop: no type 2! After time " << timeIndex << endl;
                    coexistence = false;
                    k = nbrWritting;
                    fraction = 0;
                }
                else if (popTotal.getphi(0) == phi2)
                {
                    cout << "stop: no type 1! After time " << timeIndex << endl;
                    coexistence = false;
                    k = nbrWritting;
                    fraction = 1;
                }
                else
                {
                    cout << "stop: problem in code" << timeIndex << endl;
                    abort();
                }

            }
            if (popTotal.sizePopM() > normalParamK * 10000) //population too big: problem in code
            {
                cout << "stop: too much microbes! After time " << timeIndex << endl;
                coexistence = false;
                k = nbrWritting;
            }

        }
        popTotal.writeFile2types(nameFile,timeIndex,normalParamK,weightEnzyme,weightCarbonC,weightCarbonD,phi1,phi2);
        popTotal.writeFile2typesTotalAccess(nameFileTotal,timeIndex,normalParamK,weightEnzyme,weightCarbonC,weightCarbonD,phi1,phi2);
    }

    if ((fraction == -1) && (popTotal.nbrTypes() == 2)) fraction = popTotal.sizeSubPopM(1)/popTotal.sizePopM();

    //ofstream flux2(nameFile.c_str(),ios::app);
    //flux2 << fraction << " ";
    //for (int k(0); k < 5*squareLength*squareLength; k++)
    //    flux2 << 0 << " ";
    //flux2 << endl; //endl pour passer à la ligne

    double fractionInitialtype2((double)(nbrInitialCasesMutants)/( (double)(squareLength*squareLength) ));
    double popTotalFinal((double)popTotal.sizePopM());
    double growthRatecond(-1000), growthRate(0);
    if (fraction != 0)
    {
        growthRatecond = 1/timeIndex*log(fraction*popTotalFinal/
                                                    (double)(nbrInitialCasesMutants*nbrInitialMParCase));
        growthRate = growthRatecond;
    }
    else
    {
        growthRate = 1/timeIndex*log(1/(double)(nbrInitialCasesMutants*nbrInitialMParCase));

    }

    ostringstream oss2;
    oss2 << "./communs/scoresCommuns" << argv[1] << "_" << argv[2] << ".txt";
    string nameFileCommun(oss2.str());
    ofstream fluxcommun(nameFileCommun.c_str());
    fluxcommun << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5];
    fluxcommun << " " << argv[6] << " " << argv[7] << " " << argv[8] << " " << argv[9] << " " << argv[10];
    fluxcommun << " " << argv[11] << " " << argv[12] << " " <<  argv[13] << " " <<  argv[14] << " " <<  argv[15];
    fluxcommun << " " <<  argv[16] << " " << argv[17] << " " << argv[18] << " " << fractionInitialtype2;
    fluxcommun << " " << fraction << " " << timeIndex << " " << popTotalFinal << " " << growthRate;
    fluxcommun << " " << growthRatecond << endl;

    return 0;
}
