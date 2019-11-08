#include <iostream>
#include "./code/Population.h"
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char* argv[])
{


    //Argument
    cout << "argc = " << argc << endl;
    if (argc != 8)
        {
            cout << "7 arguments: ID cluster, ID job, phi (phi1),";
            cout << " K, timeMax, probaMutation";
            cout << "maximum interval between 2 calculations of birth" << endl;
            return -1;
        }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //parameters to specify
    //initialisation of the stream
    srand(atof(argv[1])+atof(argv[2]));


    //parameter of normalisation K
    //The idea of this program is to be able to look at large K values, typically K~alpha=wM/wD, which we can't do if the 4 species are random because the simulations would be too long.
    double const normalParamK(atof(argv[4])); // -> mean number of individuals

    //parameters of the entities:
    double const erosion(0);
    //parameters of cost:
    double const weightMicrobe(1e-9);  //wM -> biomass of 1 microbe (mg)
    double const weightCarbonC(1e-16);  //wC -> biomass of 1 carbonC (mg) = typically cellulose
    double const weightEnzyme(1e-16);  //wZ -> biomass of 1 enzyme (mg) = typically cellulase
    double const weightCarbonD(1e-19);  //wD -> biomasse of 1 carbonD (mg) = typically glucose
    //we need to have alpha=wM/wD > beta=wC/wD > gamma=wZ/wD

    //Microbe
    double const VmU(0.42);
    double const KmU(3e-10);
    double const indDeathRateM(2e-4); //dM
    double const gammaM(0.3);
    double const gammaZ(0.4);
    double const paramP(0.5);  //p
    double const phi1(atof(argv[3]));  //initial value of phi

    //Enzyme
    double const indDeathRateZ(2e-3);  //dZ

    //CarbonC
    double const IC(5e-13);
    double const indDeathCRate(1e-6);  //l_C
    double const VmD(7e+5);


    //CarbonD
    double const ID(0);
    double const indDeathDRate(1e-6);  //l_D


    //mutation
    double const probaMutation(atof(argv[6]));  //probability of mutation at each event of birth.
    double const paramMutation(0.02); //standard deviation of the Gaussian r.v.


    //Time and file parameters to be adjusted
    double const nbrWritting(1000);  //Numbers of intervals time to write on the file
    double const timeMax(atof(argv[5]));  //Maximal time to observe
    double const timeStepEuler(0.001);  //step for solving differential equations
    double const annexTimeStep(atof(argv[7]));

  
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
    int nbrInitialMicrobes(0);  // Initial number of Microbes -> en nbr d'individus total par case
    double sizeInitialZ(0); //(biomasseInitialZ/(normalParamK*weightEnzyme));
    double sizeInitialC(0); //(biomasseInitialC/(normalParamK*weightCarbonC));
    double sizeInitialD(0); //(biomasseInitialD/(normalParamK*weightCarbonD));
    double sizeInitialCbis(atof(argv[5]));

    bool stop(false);



    eqAttendus(phi1, gammaM, VmU, indDeathRateM, gammaZ, indDeathRateZ, KmU, VmD, indDeathCRate, indDeathDRate,
                 erosion, IC, ID, paramP, weightMicrobe, weightEnzyme, weightCarbonC, weightCarbonD, normalParamK,
                 nbrInitialMicrobes, sizeInitialZ, sizeInitialC, sizeInitialD);


    Microbe microbeInit(structuralCostCarbonDMicrobeNormal,energeticCostCarbonDMicrobeNormal,phi1,VmU,paramH,indDeathRateM,gammaM,
                        bargammaZ,paramP);


    Population popTotal(sizeInitialZ,sizeInitialC,sizeInitialD,indDeathRateZ,structuralCostCarbonDEnzyme,
                        energeticCostCarbonDEnzyme,indDeathCRate,creationCRate,indBindCZRate,structuralCostCarbonDCarbonC,
                        indDeathDRate,creationDRate,erosion,paramH,nbrInitialMicrobes,microbeInit);


    //other time and file parameters
    double timeIndex(0);
    double timeStep(0);

    popTotal.writeFile(nameFileTotal,timeIndex,normalParamK,weightEnzyme,weightCarbonC,weightCarbonD);


    for (int k(1); k <= nbrWritting; k++)
    {
        while (timeIndex < (timeMax / nbrWritting)*k && !stop)
        {
            popTotal.oneStepAlgo(timeStep,timeStepEuler,probaMutation,paramMutation,annexTimeStep);
            timeIndex += timeStep;


            if (popTotal.sizePopM() == 0) // no more individuals
            {
                cout << "stop: plus de bacteries " << timeIndex << endl;
                stop = true;
                k = nbrWritting;

            }
            if (popTotal.sizePopM() > normalParamK * 10000) //too big population --> problem in code
            {
                cout << "stop: too much microbes! After time " << timeIndex << endl;
                stop = true;
                k = nbrWritting;
            }

        }
        popTotal.writeFile(nameFileTotal,timeIndex,normalParamK,weightEnzyme,weightCarbonC,weightCarbonD);
    }




    double popTotalFinal((double)popTotal.sizePopM());

    ostringstream oss2;
    oss2 << "./communs/scoresCommuns" << argv[1] << "_" << argv[2] << ".txt";
    string nameFileCommun(oss2.str());
    ofstream fluxcommun(nameFileCommun.c_str());
    fluxcommun << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6] << " ";
    fluxcommun << argv[7] << " ";
    fluxcommun << timeIndex << " " << popTotalFinal << " " << endl;

    return 0;
}
