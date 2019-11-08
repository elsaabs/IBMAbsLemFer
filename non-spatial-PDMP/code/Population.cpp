#include "Population.h"

using namespace std;



/*functions outside of the class*/
//uniform random variable between 0 and maxInt-1
int varUnif(int const& maxInt)
{
    return (rand() % maxInt);
}

//uniform random variable on [0,1]
double varUniform()
{
    return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

//exponential random variable with parameter "parameter"
double varExpo(double const& parameter)
{
    return -log( 1.0 - ((double)rand()+1.0)/((double)RAND_MAX+2.0)) / parameter;
}


void eqAttendus(double const& phi, double const& gammaM, double const& VmU, double const& dM,
                  double const& gammaZ, double const& dZ, double const& KmU, double const& VmD, double const& lC, double const& lD,
                  double const& erosion, double const& IC, double const& ID, double const& paramP, double const& weightMicrobe,
                  double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD, double const& normalParamK,
                  int &MeqInd, double &ZeqDens, double &CeqDens, double &DeqDens)
{
    if ((1-phi)*gammaM*VmU-dM < 0)
    //No steady state with M positive
    {
        cout << "No steady state with M positive"; abort();
    }
    else
    {
        double C1( phi*gammaZ*dM / ((1-phi)*gammaM*dZ) );
        double deq( KmU*dM / ((1-phi)*gammaM*VmU-dM) );
        double C2( (1-erosion)*((1-paramP)*dM + dZ*C1) - dM / ((1-phi)*gammaM) );
        double Apoly( VmD*C1*( (1-erosion)*paramP*dM + C2 ) );
        double Bpoly( VmD*C1*( IC+ID-lD*deq ) + lC*C2 );
        double Cpoly( ( ID - lD*deq )*lC );
        double Disc( Bpoly*Bpoly - 4*Apoly*Cpoly );
        double coef( 1-(1-erosion)*(phi*gammaZ-(1-phi)*gammaM) );
        double ceq(0);
        double sgndeterminant(0);
        if (Disc < 0) {cout << "No steady state with M positive"; abort();}
        else
        {
            double m1( (-Bpoly-sqrt(Disc)) / (2*Apoly) );
            double m2( (-Bpoly+sqrt(Disc)) / (2*Apoly) );
            double meq( max(m1,m2) );
            if ( meq < 0 ) {cout << "No steady state with M positive"; abort();}
            else
            {
            meq = m1;
            ceq = ((IC+(1-erosion)*paramP*dM*meq) / (lC+VmD*C1*meq));
            sgndeterminant = (lC+VmD*C1*meq)*coef*dM / ((1-phi)*gammaM) + (1-erosion)*paramP*dM*lC - C1*lC*VmD*ceq;
            if ( sgndeterminant >= 0 )
            {
            MeqInd = max(floor(meq*normalParamK/weightMicrobe),1.0);
            ZeqDens = C1*meq/weightEnzyme;
            CeqDens = ceq/weightCarbonC;
            DeqDens = deq/weightCarbonD;
            //cout << "sign of the determinant ~ " << sgndeterminant << endl << endl;
            }
            else
            {
            meq = m2;
            ceq = ((IC+(1-erosion)*paramP*dM*meq) / (lC+VmD*C1*meq));
            sgndeterminant = (lC+VmD*C1*meq)*coef*dM / ((1-phi)*gammaM) + (1-erosion)*paramP*dM*lC - C1*lC*VmD*ceq;
            if ( sgndeterminant >= 0 )
            {
                MeqInd = max(floor(meq*normalParamK/weightMicrobe),1.0);
                ZeqDens = C1*meq/weightEnzyme;
                CeqDens = ceq/weightCarbonC;
                DeqDens = deq/weightCarbonD;
            }
            else {cout << "No steady state with M positive"; abort();}
            }
            }
        }
        }
}






/*---------------------------------------------------------*/

/*functions of the class: Population*/




Population::Population(double const& sizeInitialEnzyme, double const& sizeInitialC,
                 double const& sizeInitialD, double const& indDeathZRate, double const& nbrEquivZD, double const& energCostZD,
                 double const& indDeathCRate, double const& creationCRate, double const& indBindZCRate,
                 double const& nbrEquivCD, double const& indDeathDRate, double const& creationDRate,
                 double const& erosion, double const& paramH, int const& nbrInitialMicrobes, Microbe const& microbeInitial):
                     m_indDeathZRate(indDeathZRate), m_nbrEquivZD(nbrEquivZD), m_energCostZD(energCostZD),
                     m_indDeathCRate(indDeathCRate), m_creationCRate(creationCRate), m_indBindZCRate(indBindZCRate),
                     m_nbrEquivCD(nbrEquivCD), m_indDeathDRate(indDeathDRate), m_creationDRate(creationDRate),
                     m_erosion(erosion), m_paramH(paramH) 
{
    m_linearGrowthtempInd = 0;


    //initialization of types Z,C,D
    m_sizeZ = sizeInitialEnzyme;
    m_sizeC = sizeInitialC;
    m_sizeD = sizeInitialD;


    //initialization of bacterial individuals
    vector<double> newTabPop(nbrInitialMicrobes,0);
    for (int i(0); i < nbrInitialMicrobes; i++)
    {
        newTabPop[i] = varUniform();
    }
    vector<vector<double> > newInd(1,newTabPop);
    m_individuals = newInd;

    vector<Microbe> newTabMicrobe(1,microbeInitial);
    m_microbeTypes = newTabMicrobe;
    m_sizeMTotal = nbrInitialMicrobes;


    m_productionZ = (double)nbrInitialMicrobes*microbeInitial.maxproduceZRate();
    m_maxlinearGrowth = (double)nbrInitialMicrobes*microbeInitial.maxlineargrowthM();
    m_gammaZ = microbeInitial.getgammaZ();
    m_maxCostD = microbeInitial.getmaxCost();
    m_maxGrowthrateforanyphi = microbeInitial.getmaxGrowthforanyphi();


    m_maxdeathrate = microbeInitial.maxdeathRate();
}








Population::~Population()
{
    //dtor
}




//related to Z,C,D

void Population::oneStepEuler(double const& sizeStep)
{
    m_linearGrowthtempInd = m_linearGrowthtempInd
                                   + m_maxGrowthrateforanyphi * ( m_sizeD / ( m_paramH + m_sizeD ) ) * sizeStep;

    m_sizeZ = m_sizeZ + (m_productionZ * ( m_sizeD / ( m_paramH + m_sizeD ))
                         - m_indDeathZRate * m_sizeZ) * sizeStep;

    m_sizeC = m_sizeC + (m_creationCRate - m_indDeathCRate * m_sizeC
                         - m_indBindZCRate * m_sizeC * m_sizeZ) * sizeStep;

    m_sizeD = m_sizeD + (m_creationDRate - m_indDeathDRate * m_sizeD
                         + m_indBindZCRate * m_sizeC * m_sizeZ * m_nbrEquivCD
                         + (1-m_erosion) * m_indDeathZRate * m_sizeZ * m_nbrEquivZD
                         - (m_productionZ/m_gammaZ + m_maxlinearGrowth * m_maxCostD)
                                       * ( m_sizeD / ( m_paramH + m_sizeD ))) * sizeStep;
}




void Population::eulerScheme(double const& sizeStep, double const& intervalTime, vector<double> &birthTab)
{
    birthTab.clear();
    m_linearGrowthtempInd = 0;
    
    double time(0);
    while (time < intervalTime)
    {
        if (time+sizeStep <= intervalTime)
        {
            oneStepEuler(sizeStep);
            time += sizeStep;
        }
        else
        {
            oneStepEuler(intervalTime-time);
            time = intervalTime;
        }
    }

    for (int k(0); k < m_individuals.size(); k++)
    {
        double phitemp(m_microbeTypes[k].getphi());
        for (int k2(0); k2 < m_individuals[k].size(); k2++)
        {
            m_individuals[k][k2] += (1-phitemp)*m_linearGrowthtempInd;
            double growth(m_individuals[k][k2]);
            if (growth > 1)
            {
                if (growth > 1.01) cout << "Time steps too large " << endl;
                birthTab.push_back((growth-1)/2);
                birthTab.push_back(k);
                birthTab.push_back(k2);
                m_individuals[k][k2] = (growth-1)/2;
            }
        }
    }

}







double Population::totalZ() const
{
    return m_sizeZ;
}


double Population::totalC() const
{
    return m_sizeC;
}

double Population::totalD() const
{
    return m_sizeD;
}




//related to M


int Population::nbrTypes() const
{
    return m_individuals.size();
}




double Population::sizePopM() const
{
    return m_sizeMTotal;
}





double Population::sizeSubPopM(int const& indexType) const
{
    if (indexType >= m_individuals.size())
    {
        cout << "problem with function sizeSubPop" << endl;
        assert(false);
    }
    else
    {
        return (double)m_individuals[indexType].size();
    }
}





void Population::functionIndexM(double const& indexMicrobe, int &indexSubPopM, int &indexInPopM) const
{
    if (indexMicrobe < 0) cout << "Problem: no microbe..." << endl;
    indexInPopM = indexMicrobe;
    indexSubPopM = 0;
    bool test(true);
    while (test)
    {
        if (indexInPopM < sizeSubPopM(indexSubPopM))
        {
            test = false;
        }
        else
        {
            indexInPopM -= sizeSubPopM(indexSubPopM);
            indexSubPopM += 1;
        }
    }
}










double Population::maxDeathRate()
{
    double typeRate(0);
    double maxRate(0);
    if (m_maxdeathrate == -1)
    {
        for (int j(0); j<nbrTypes(); j++)
        {
            typeRate = m_microbeTypes[j].maxdeathRate();
            if (typeRate > maxRate) maxRate = typeRate;
        }
        m_maxdeathrate = maxRate;
    }
    if (m_maxdeathrate < -0.00001) {cout << "problem with maxdeathrate"; abort();}
    return m_maxdeathrate;
}






double Population::deathRate(double const& indexMicrobe) const
{
    int indexSubPopM(0), indexInPopM(0);
    functionIndexM(indexMicrobe,indexSubPopM,indexInPopM);
    return m_microbeTypes[indexSubPopM].deathRate();
}










void Population::addMaxRate(int const& indexSubPopM)
{
    if (m_microbeTypes[indexSubPopM].maxdeathRate() > m_maxdeathrate)
            m_maxdeathrate = m_microbeTypes[indexSubPopM].maxdeathRate();
}



void Population::suppMaxRate(double const& death)
{
    if (death >= m_maxdeathrate-0.000001)
        {
            m_maxdeathrate = -1;
            maxDeathRate();
        }
}






void Population::updateTotals(int const& indexSubPopM, int const& indexInPopM, bool const& isplus)
{
    if (isplus)
    {
        m_maxlinearGrowth += m_microbeTypes[indexSubPopM].maxlineargrowthM();
        m_productionZ += m_microbeTypes[indexSubPopM].maxproduceZRate();
    }
    else
    {
        m_maxlinearGrowth -= m_microbeTypes[indexSubPopM].maxlineargrowthM();
        m_productionZ -= m_microbeTypes[indexSubPopM].maxproduceZRate();
    }
}








void Population::addMicrobe(double const& initialGrowth, int const& indexSubPopM, int const& indexInPopM,
                            double const& probaMutation, double const& paramMutation)
{
    if (varUniform() < probaMutation)
    {
        double oldPhi(m_microbeTypes[indexSubPopM].getphi());
        double newPhi(0);
        do
            {
                newPhi = oldPhi + paramMutation * sqrt(-2.0 * log( ((double)rand()+1.0)/((double)RAND_MAX+2.0) ))
                         * cos( 6.2852 * (double)rand() /((double)RAND_MAX+1.0) );
            }
        while ((newPhi<=0)||(newPhi>=1));
        Microbe newMicrobe(m_microbeTypes[indexSubPopM],newPhi);
        m_microbeTypes.push_back(newMicrobe);
        vector<double> newvector(1,initialGrowth);
        m_individuals.push_back(newvector);
        addMaxRate(m_microbeTypes.size()-1);
        updateTotals(m_microbeTypes.size()-1,0, true);
    }
    else
    {
        m_individuals[indexSubPopM].push_back(initialGrowth);
        updateTotals(indexSubPopM,m_individuals[indexSubPopM].size()-1,true);
    }
    m_sizeMTotal += 1;

}













void Population::deleteMicrobe(double &carbonCReleased, double &carbonDReleased, double const& indexMicrobe)
{
    int indexSubPop(0), indexInPop(0);
    if (indexMicrobe < 0) cout << "problem when deleting microbe: no microbe!" << endl;
    functionIndexM(indexMicrobe,indexSubPop,indexInPop);
    double actualGrowth(m_individuals[indexSubPop][indexInPop]);
    m_microbeTypes[indexSubPop].carbonReleasedwhenDeath(carbonCReleased, carbonDReleased, m_nbrEquivCD, actualGrowth);
    if (sizeSubPopM(indexSubPop) == 1)
    {
        updateTotals(indexSubPop,indexInPop, false);
        double death(m_microbeTypes[indexSubPop].maxdeathRate());
        m_individuals.erase(m_individuals.begin()+indexSubPop);
        m_microbeTypes.erase(m_microbeTypes.begin()+indexSubPop);
        suppMaxRate(death);
    }
    else
    {
        updateTotals(indexSubPop,indexInPop, false);
        m_individuals[indexSubPop].erase(m_individuals[indexSubPop].begin()+indexInPop);
    }
    m_sizeMTotal -= 1;
}







void Population::oneStepAlgo(double &timeStep, double const& timeStepEuler, double const& probaMutation,
                                    double const& paramMutation, double const& annexTimeStep)
{
    double microbeDeathRate(maxDeathRate());
    double va1(varUniform() * microbeDeathRate);
    if (m_sizeMTotal == 0) {cout << "probleme : oneStep alors que pas de bacteries..."; abort();}
    timeStep = varExpo(microbeDeathRate*m_sizeMTotal); //gives the time before the next death
    if (timeStep < 0) {cout << "probleme avec le calcul de timeStep"; abort();}


    vector<double> birthTab;

    if (timeStep >= annexTimeStep)
    {
        timeStep = annexTimeStep;
        eulerScheme(timeStepEuler, timeStep, birthTab);
        for (int k(0); k<birthTab.size()/3; k++)
            addMicrobe(birthTab[3*k],(int)birthTab[3*k+1],(int)birthTab[3*k+2],probaMutation,paramMutation);
    }
    else
    {
        int vaIndex(varUnif(m_sizeMTotal));
        int indexSubPop(0), indexInPop(0);
        functionIndexM(vaIndex,indexSubPop,indexInPop);

        eulerScheme(timeStepEuler, timeStep, birthTab);
        for (int k(0); k<birthTab.size()/3; k++)
            addMicrobe(birthTab[3*k],(int)birthTab[3*k+1],(int)birthTab[3*k+2],probaMutation,paramMutation);
       
        if (va1 < m_microbeTypes[indexSubPop].deathRate())
        {
            double carbonCReleased(0), carbonDReleased(0);
            deleteMicrobe(carbonCReleased,carbonDReleased,vaIndex);
            m_sizeC += (1-m_erosion)*carbonCReleased;
            m_sizeD += (1-m_erosion)*carbonDReleased;
        }
    }


}


double Population::getphi(int const& indexSubPop) const
{
    return m_microbeTypes[indexSubPop].getphi();
}





void Population::writeFile(std::string nameFile, double const& timeIndex, double const& normalParamK,
                           double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD)
{
    ofstream flux(nameFile.c_str(),ios::app);
    for (int k(0); k < m_microbeTypes.size(); k++)
    {
    flux << timeIndex << " ";
    flux << m_sizeZ*normalParamK*weightEnzyme << " ";
    flux << m_sizeC*normalParamK*weightCarbonC << " ";
    flux << m_sizeD*normalParamK*weightCarbonD << " ";
    flux << m_sizeMTotal << " " << m_microbeTypes.size() << " ";
    flux << m_microbeTypes[k].getphi() << " ";
    flux << m_individuals[k].size();
    flux << endl;
    }
}






void Population::writeFile2types(std::string nameFile, double const& timeIndex, double const& normalParamK,
                           double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD,
                           double const& phi1, double const& phi2)
{
    ofstream flux(nameFile.c_str(),ios::app);
    flux << timeIndex << " ";
    flux << m_sizeZ*normalParamK*weightEnzyme << " ";
    flux << m_sizeC*normalParamK*weightCarbonC << " ";
    flux << m_sizeD*normalParamK*weightCarbonD << " ";
    flux << m_individuals[0].size() << " ";
    flux << m_individuals[1].size();
    flux << endl;
}










