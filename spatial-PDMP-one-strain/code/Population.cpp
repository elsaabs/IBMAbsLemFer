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
                    else {cout << "Pas d'equilibre avec M positif."; abort();}
                }
            }
        }
    }
}






/*---------------------------------------------------------*/

/*functions de la classe Population*/






// construction for 1 type/strain of individuals that fill up the whole lattice
Population::Population(int const& squareLength, double const& sizeInitialEnzyme, double const& sizeInitialC,
                       double const& sizeInitialD,
                       double const& indDeathZRate, double const& nbrEquivZD, double const& energCostZD,
                       double const& indDeathCRate, double const& creationCRate, double const& indBindZCRate,
                       double const& nbrEquivCD, double const& indDeathDRate, double const& creationDRate,
                       double const& erosion, double const& paramH, double const& diffusionD,
                       int const& nbrParCase, Microbe const& microbeInitial):
m_length(squareLength), m_indDeathZRate(indDeathZRate), m_nbrEquivZD(nbrEquivZD), m_energCostZD(energCostZD),
m_indDeathCRate(indDeathCRate), m_creationCRate(creationCRate), m_indBindZCRate(indBindZCRate),
m_nbrEquivCD(nbrEquivCD), m_indDeathDRate(indDeathDRate), m_creationDRate(creationDRate),
m_erosion(erosion), m_paramH(paramH), m_diffusionD(diffusionD)
{
    vector<double> newtabGrowth(squareLength*squareLength,0);
    m_linearGrowthtempInd = newtabGrowth;

    //initialization of types Z,C,D
    vector<double> newtabZ(squareLength*squareLength,sizeInitialEnzyme);
    m_sizeZ = newtabZ;

    vector<double> newtabC(squareLength*squareLength,sizeInitialC);
    m_sizeC = newtabC;

    vector<double> newtabD(squareLength*squareLength,sizeInitialD);
    m_sizeD = newtabD;


    //initialization of bacteria
    vector<double> newTabPop(nbrParCase*squareLength*squareLength,0);
    vector<int> newTabPosition(nbrParCase*squareLength*squareLength,0);
    for (int i(0); i < nbrParCase*squareLength*squareLength; i++)
    {
        newTabPosition[i] = (int)(i/nbrParCase);
        newTabPop[i] = varUniform();
    }
    vector<vector<double> > newInd(1,newTabPop);
    m_individuals = newInd;
    vector<vector<int> > newPosition(1,newTabPosition);
    m_positions = newPosition;

    vector<int> newTabsizeM(squareLength*squareLength,nbrParCase);
    m_sizeM = newTabsizeM;

    vector<Microbe> newTabMicrobe(1,microbeInitial);
    m_microbeTypes = newTabMicrobe;
    m_sizeMTotal = nbrParCase*squareLength*squareLength;

    vector<double> total(squareLength*squareLength,(double)nbrParCase*microbeInitial.maxproduceZRate());
    vector<double> total2(squareLength*squareLength,(double)nbrParCase*microbeInitial.maxlineargrowthM());
    m_productionZ = total;
    m_maxlinearGrowth = total2;
    m_gammaZ = microbeInitial.getgammaZ();
    m_maxCostD = microbeInitial.getmaxCost();
    m_maxGrowthrateforanyphi = microbeInitial.getmaxGrowthforanyphi();

    //m_maxabsorptionrate = microbeInitial.maxabsorptionRate();
    m_maxdeathrate = microbeInitial.maxdeathRate();


    vector<double> types(squareLength*squareLength,microbeInitial.getphi());
    m_types = types;

}






//3rd construction: initialisation -> microsites with mutants at the mutant phi value steady state, microsites with residents at the resident phi value steady state

void Population::change(double const& sizeInitialEnzyme, double const& sizeInitialC, double const& sizeInitialD,
                        double const& sizeInitialEnzyme2, double const& sizeInitialC2, double const& sizeInitialD2,
                        int const& nbrParCase, int const& nbrParCase2, int const& numeroCaseInitiale, int const& nbrCasesMutant,
                        Microbe const& microbe1Initial, Microbe const& microbe2Initial)
{
    //check number of microsites with mutants
    if (nbrCasesMutant > m_length*m_length) {cout << "On se concentre: trop de cases de mutants..." << endl; abort();}

    //construction indexes of microsites with mutants
    int j(0), k(1), kbis(1), indiceSigne(1), indiceAdd(1);
    if (numeroCaseInitiale >= m_length*m_length) {cout << "NumeroCaseInitale trop grand" << endl; abort();}
    int casesMutant[nbrCasesMutant];
    casesMutant[0] = numeroCaseInitiale;
    for (int l(1); l < nbrCasesMutant; l++)
    {
        if (indiceSigne == 1) casesMutant[l] = casesMutant[l-1] + indiceAdd;
        if (indiceSigne == 0) casesMutant[l] = casesMutant[l-1] - indiceAdd;
        if (l == k)
        {
            if (indiceAdd == 1) {indiceAdd = m_length; k += kbis;}
            else {indiceAdd = 1; kbis++; k += kbis; indiceSigne = (indiceSigne+1)%2;}
        }
        if (casesMutant[l] < 0 || casesMutant[l] >= m_length*m_length) {cout << "probleme calcul casesMutant" << endl; abort();}
    }

    //construction indexes of microsites with residents
    int casesResidents[m_length*m_length-nbrCasesMutant];
    bool test(true);
    k = 0;
    for (int l(0); l < m_length*m_length; l++)
    {
        test = true;
        for (int ll(0); ll < nbrCasesMutant; ll++)
        {
            if (l == casesMutant[ll]) test = false;
        }
        if (test) {casesResidents[k] = l; k++;}
    }
    if (k > m_length*m_length-nbrCasesMutant)  {cout << "probleme calcul casesResidents" << endl; abort();}


    //initialization carbon
    vector<double> newtabZ(m_length*m_length,sizeInitialEnzyme);
    for (int l(0); l < nbrCasesMutant; l++)
    {
        newtabZ[casesMutant[l]] = sizeInitialEnzyme2;
    }
    m_sizeZ = newtabZ;

    vector<double> newtabC(m_length*m_length,sizeInitialC);
    for (int l(0); l < nbrCasesMutant; l++)
    {
        newtabC[casesMutant[l]] = sizeInitialC2;
    }
    m_sizeC = newtabC;

    vector<double> newtabD(m_length*m_length,sizeInitialD);
    for (int l(0); l < nbrCasesMutant; l++)
    {
        newtabD[casesMutant[l]] = sizeInitialD2;
    }
    m_sizeD = newtabD;


    //initialization of bacteria
    vector<double> newTabPop1(nbrParCase*(m_length*m_length-nbrCasesMutant),0);
    vector<int> newTabPosition1(nbrParCase*(m_length*m_length-nbrCasesMutant),0);
    for (int i(0); i < nbrParCase*(m_length*m_length-nbrCasesMutant); i++)
    {
        newTabPop1[i] = varUniform();
        newTabPosition1[i] = casesResidents[i/nbrParCase];
    }
    vector<double> newTabPop2(nbrParCase2*nbrCasesMutant,0);
    vector<int> newTabPosition2(nbrParCase2*nbrCasesMutant,0);
    for (int i(0); i < nbrParCase2*nbrCasesMutant; i++)
    {
        newTabPosition2[i] = casesMutant[i/nbrParCase2];
        newTabPop2[i] = varUniform();
    }
    vector<vector<double> > newInd(1,newTabPop1);
    newInd.push_back(newTabPop2);
    m_individuals = newInd;
    vector<vector<int> > newPosition(1,newTabPosition1);
    newPosition.push_back(newTabPosition2);
    m_positions = newPosition;

    vector<Microbe> newTabMicrobe(1,microbe1Initial);
    newTabMicrobe.push_back(microbe2Initial);
    m_microbeTypes = newTabMicrobe;
    m_sizeMTotal = nbrParCase*(m_length*m_length-nbrCasesMutant)+nbrParCase2*nbrCasesMutant;


    vector<double> total(m_length*m_length,(double)nbrParCase*microbe1Initial.maxproduceZRate());
    for (int l(0); l < nbrCasesMutant; l++)
    {
        total[casesMutant[l]] = (double)nbrParCase2*microbe2Initial.maxproduceZRate();
    }

    vector<double> total2(m_length*m_length,(double)nbrParCase*microbe1Initial.maxlineargrowthM());
    for (int l(0); l < nbrCasesMutant; l++)
    {
        total2[casesMutant[l]] = (double)nbrParCase2*microbe2Initial.maxlineargrowthM();
    }

    m_productionZ = total;
    m_maxlinearGrowth = total2;
    m_gammaZ = microbe1Initial.getgammaZ();
    m_maxCostD = microbe1Initial.getmaxCost();
    m_maxGrowthrateforanyphi = microbe1Initial.getmaxGrowthforanyphi();


    m_maxdeathrate = max(microbe1Initial.maxdeathRate(),microbe2Initial.maxdeathRate());



    vector<double> types(m_length*m_length,microbe1Initial.getphi());
    for (int l(0); l < nbrCasesMutant; l++)
    {
        types[casesMutant[l]] = microbe2Initial.getphi();
    }
    m_types = types;

}





Population::~Population()
{
    //dtor
}



//related to Z,C,D

void Population::oneStepEulerOneCase(double const& sizeStep, int const& index)
{
    m_linearGrowthtempInd[index] = m_linearGrowthtempInd[index]
    + m_maxGrowthrateforanyphi * ( m_sizeD[index] / ( m_paramH + m_sizeD[index] ) ) * sizeStep;

    m_sizeZ[index] = m_sizeZ[index] + (m_productionZ[index] * ( m_sizeD[index] / ( m_paramH + m_sizeD[index] ))
                                       - m_indDeathZRate * m_sizeZ[index]) * sizeStep;

    m_sizeC[index] = m_sizeC[index] + (m_creationCRate - m_indDeathCRate * m_sizeC[index]
                                       - m_indBindZCRate * m_sizeC[index] * m_sizeZ[index]) * sizeStep;

    m_sizeD[index] = m_sizeD[index] + (m_creationDRate - m_indDeathDRate * m_sizeD[index]
                                       + m_indBindZCRate * m_sizeC[index] * m_sizeZ[index] * m_nbrEquivCD
                                       + (1-m_erosion) * m_indDeathZRate * m_sizeZ[index] * m_nbrEquivZD
                                       - (m_productionZ[index]/m_gammaZ + m_maxlinearGrowth[index] * m_maxCostD)
                                       * ( m_sizeD[index] / ( m_paramH + m_sizeD[index] ))) * sizeStep;
}



void Population::eulerSchemeOneCase(double const& sizeStep, double const& intervalTime, int const& index)
{
    if (index > m_length*m_length)
    {
        cout << "problem with eulerSchemeOneCase";
        abort();
    }
    else
    {
        double time(0);
        while (time < intervalTime)
        {
            if (time+sizeStep <= intervalTime)
            {
                oneStepEulerOneCase(sizeStep, index);
                time += sizeStep;
            }
            else
            {
                oneStepEulerOneCase(intervalTime-time, index);
                time = intervalTime;
            }
        }
    }
}



void Population::eulerScheme(double const& sizeStep, double const& intervalTime, vector<double> &birthTab)
{
    vector<double> newTabZero(m_length*m_length,0);
    birthTab.clear();
    m_linearGrowthtempInd = newTabZero;
    for (int i(0); i<m_length*m_length; i++)
    {
        eulerSchemeOneCase(sizeStep,intervalTime,i);
    }

    for (int k(0); k < m_individuals.size(); k++)
    {
        double phitemp(m_microbeTypes[k].getphi());
        for (int k2(0); k2 < m_individuals[k].size(); k2++)
        {
            m_individuals[k][k2] += (1-phitemp)*m_linearGrowthtempInd[m_positions[k][k2]];
            double growth(m_individuals[k][k2]);
            if (growth > 1)
            {
                if (growth > 1.01) cout << "pasTropGd " << endl;
                birthTab.push_back((growth-1)/2);
                birthTab.push_back(k);
                birthTab.push_back(m_positions[k][k2]);
                m_individuals[k][k2] = (growth-1)/2;
            }
        }
    }

}









void Population::diffusionCarbon(double const& sizeTimeStep)
{
    double coef(sizeTimeStep*m_diffusionD);
    if (coef > 0.2) {coef = 0.20;}
    int j(0),k(0);
    vector<double> sizeDtemp(m_length*m_length,0);
    sizeDtemp = m_sizeD;
    //top-left
    sizeDtemp[0] += coef * (m_sizeD[1] + m_sizeD[m_length] - 2*m_sizeD[0]);
    //top-right
    sizeDtemp[m_length - 1] += coef * (m_sizeD[m_length-2] + m_sizeD[2*m_length-1] - 2*m_sizeD[m_length-1]);
    //bottom-left
    sizeDtemp[m_length*(m_length-1)] += coef * (m_sizeD[m_length*(m_length-2)] + m_sizeD[m_length*(m_length-1)+1]
                                                - 2*m_sizeD[m_length*(m_length-1)]);
    //bottom-right
    sizeDtemp[m_length*m_length-1] += coef * (m_sizeD[m_length*(m_length-1)-1] + m_sizeD[m_length*m_length-2]
                                              - 2*m_sizeD[m_length*m_length-1]);
    //if on an eged of the lattice
    //top edge
    for (k = 1; k<m_length-1; k++)
    {sizeDtemp[k] += coef * (m_sizeD[k+1] + m_sizeD[k-1] + m_sizeD[k+m_length] - 3*m_sizeD[k]);}
    //right edge
    for (j = 1; j<m_length-1; j++)
    {k = (j+1)*m_length - 1;
        sizeDtemp[k] += coef * (m_sizeD[k-m_length] + m_sizeD[k-1] + m_sizeD[k+m_length] - 3*m_sizeD[k]);}
    //bottom edge
    for (k = m_length*(m_length-1)+1; k<m_length*m_length-1; k++)
    {sizeDtemp[k] += coef * (m_sizeD[k-1] + m_sizeD[k+1] + m_sizeD[k-m_length] - 3*m_sizeD[k]);}
    //left edge
    for (j = 1; j<m_length-1; j++)
    {k = j*m_length;
        sizeDtemp[k] += coef * (m_sizeD[k-m_length] + m_sizeD[k+1] + m_sizeD[k+m_length] - 3*m_sizeD[k]);}
    //if inside
    for (int i(1); i<m_length-1; i++)
    {
        for (j=1; j<m_length-1; j++)
        {k = j*m_length + i;
            sizeDtemp[k] += coef * (m_sizeD[k-m_length] + m_sizeD[k+1] + m_sizeD[k+m_length] + m_sizeD[k-1] - 4*m_sizeD[k]);}
    }

    m_sizeD = sizeDtemp;
}



double Population::totalZ() const
{
    double sum(0);
    for (int k(0); k < m_length*m_length; k++) sum += m_sizeZ[k];
    return sum;
}


double Population::totalC() const
{
    double sum(0);
    for (int k(0); k < m_length*m_length; k++) sum += m_sizeC[k];
    return sum;
}

double Population::totalD() const
{
    double sum(0);
    for (int k(0); k < m_length*m_length; k++) sum += m_sizeD[k];
    return sum;
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
        m_maxlinearGrowth[m_positions[indexSubPopM][indexInPopM]] += m_microbeTypes[indexSubPopM].maxlineargrowthM();
        m_productionZ[m_positions[indexSubPopM][indexInPopM]] += m_microbeTypes[indexSubPopM].maxproduceZRate();
    }
    else
    {
        m_maxlinearGrowth[m_positions[indexSubPopM][indexInPopM]] -= m_microbeTypes[indexSubPopM].maxlineargrowthM();
        m_productionZ[m_positions[indexSubPopM][indexInPopM]] -= m_microbeTypes[indexSubPopM].maxproduceZRate();
    }
}





int Population::ifMovingOneCase(int const& indexInitial, double const& phi, double const& probaReplace)
{
    int newposition(indexInitial);
    vector<int> emptyCases;
    vector<int> fullCases;

    for (int k(0); k<4; k++)
    {
        if (k == 0) newposition = indexInitial + 1;
        if (k == 1) newposition = indexInitial + m_length;
        if (k == 2) newposition = indexInitial - 1;
        if (k == 3) newposition = indexInitial - m_length;

        if ((newposition >= 0)&&(newposition < m_length*m_length))
        {
            if (m_types[newposition] == 0) emptyCases.push_back(newposition);
            else  fullCases.push_back(newposition);
        }
    }

    if (emptyCases.size() != 0) {return emptyCases[varUnif(emptyCases.size())];}
    else
    {
        if (fullCases.size() == 0) {cout << "problem with moving" << endl; abort();}
        if (varUniform() < probaReplace)
        {
            newposition = fullCases[varUnif(fullCases.size())];
            double carbonCReleased(0);
            double carbonDReleased(0);

            //on vide la case
            vector<int> deathTab;
            deathTab.clear();
            for (int i(0); i<m_positions.size(); i++)
            {
                for (int i2(0); i2<m_positions[i].size(); i2++)
                {
                    if (m_positions[i][i2] == newposition)
                    {
                        deathTab.push_back(i);
                        deathTab.push_back(i2);
                    }
                }
            }

            for (int k((deathTab.size()-1)/2); k>=0; k--)
            {
                deleteMicrobe(carbonCReleased,carbonDReleased,deathTab[2*k],deathTab[2*k+1]);
                m_sizeC[newposition] += (1-m_erosion)*carbonCReleased;
                m_sizeD[newposition] += (1-m_erosion)*carbonDReleased;
            }
            //fin de la vidange de case


            m_types[newposition] = 0;
            return newposition;

        }
        else {return indexInitial;}
    }
}







void Population::addMicrobe(double const& initialGrowth, int const& indexSubPopM, int const& position,
                            double const& probaMoving, double const& probaReplace)
{
    int newposition(position);

    if ((m_types[position] == getphi(indexSubPopM)) || (m_types[position] == 0))
    {

        if (varUniform() < probaMoving) newposition = ifMovingOneCase(newposition, getphi(indexSubPopM), probaReplace);

        m_individuals[indexSubPopM].push_back(initialGrowth);
        m_positions[indexSubPopM].push_back(newposition);
        updateTotals(indexSubPopM,m_individuals[indexSubPopM].size()-1,true);
        m_types[newposition] = getphi(indexSubPopM);

        m_sizeMTotal += 1;
        m_sizeM[newposition] += 1;
    }
}













void Population::deleteMicrobe(double &carbonCReleased, double &carbonDReleased, double const& indexMicrobe)
{
    int indexSubPop(0), indexInPop(0);
    if (indexMicrobe < 0) cout << "problem when deleting microbe: no microbe!" << endl;
    functionIndexM(indexMicrobe,indexSubPop,indexInPop);
    double actualGrowth(m_individuals[indexSubPop][indexInPop]);
    double position(m_positions[indexSubPop][indexInPop]);
    m_microbeTypes[indexSubPop].carbonReleasedwhenDeath(carbonCReleased, carbonDReleased, m_nbrEquivCD, actualGrowth);
    if (sizeSubPopM(indexSubPop) == 1)
    {
        updateTotals(indexSubPop,indexInPop, false);
        double death(m_microbeTypes[indexSubPop].maxdeathRate());
        m_individuals.erase(m_individuals.begin()+indexSubPop);
        m_microbeTypes.erase(m_microbeTypes.begin()+indexSubPop);
        m_positions.erase(m_positions.begin()+indexSubPop);
        suppMaxRate(death);
    }
    else
    {
        updateTotals(indexSubPop,indexInPop, false);
        m_individuals[indexSubPop].erase(m_individuals[indexSubPop].begin()+indexInPop);
        m_positions[indexSubPop].erase(m_positions[indexSubPop].begin()+indexInPop);
    }
    m_sizeMTotal -= 1;
    m_sizeM[position] -= 1;
}



void Population::deleteMicrobe(double &carbonCReleased, double &carbonDReleased, int const& indexSubPop,
                               int const& indexInPop)
{
    double actualGrowth(m_individuals[indexSubPop][indexInPop]);
    double position(m_positions[indexSubPop][indexInPop]);
    m_microbeTypes[indexSubPop].carbonReleasedwhenDeath(carbonCReleased, carbonDReleased, m_nbrEquivCD, actualGrowth);
    if (sizeSubPopM(indexSubPop) == 1)
    {
        updateTotals(indexSubPop,indexInPop, false);
        double death(m_microbeTypes[indexSubPop].maxdeathRate());
        m_individuals.erase(m_individuals.begin()+indexSubPop);
        m_microbeTypes.erase(m_microbeTypes.begin()+indexSubPop);
        m_positions.erase(m_positions.begin()+indexSubPop);
        suppMaxRate(death);
    }
    else
    {
        updateTotals(indexSubPop,indexInPop, false);
        m_individuals[indexSubPop].erase(m_individuals[indexSubPop].begin()+indexInPop);
        m_positions[indexSubPop].erase(m_positions[indexSubPop].begin()+indexInPop);
    }
    m_sizeMTotal -= 1;
    m_sizeM[position] -= 1;
}








void Population::oneStepDiffusion(double &timeStep, double const& timeStepEuler, double const& probaMoving,
                                  double const& diffusionTimeStep, double const& probaReplace)
{
    double microbeDeathRate(maxDeathRate());
    double va1(varUniform() * microbeDeathRate);
    if (m_sizeMTotal == 0) {cout << "probleme : oneStep alors que pas de bacteries..."; abort();}
    timeStep = varExpo(microbeDeathRate*m_sizeMTotal); //gives the time before the next death
    if (timeStep < 0) {cout << "probleme avec le calcul de timeStep"; abort();}


    vector<double> birthTab;

    if (timeStep >= diffusionTimeStep)
    {
        timeStep = diffusionTimeStep;
        eulerScheme(timeStepEuler, timeStep, birthTab);
        for (int k(0); k<birthTab.size()/3; k++)
            addMicrobe(birthTab[3*k],(int)birthTab[3*k+1],(int)birthTab[3*k+2],probaMoving,probaReplace);
        // attention : l'enchainement de fonctions precedents creerait un probleme si il n'y avait pas exactement 2 types
        diffusionCarbon(timeStep);
    }
    else
    {
        eulerScheme(timeStepEuler, timeStep, birthTab);
        for (int k(0); k<birthTab.size()/3; k++)
            addMicrobe(birthTab[3*k],(int)birthTab[3*k+1],(int)birthTab[3*k+2],probaMoving,probaReplace);
        diffusionCarbon(timeStep);

        int vaIndex(varUnif(m_sizeMTotal));
        int indexSubPop(0), indexInPop(0);
        functionIndexM(vaIndex,indexSubPop,indexInPop);
        int positionM(m_positions[indexSubPop][indexInPop]);

        if (va1 < m_microbeTypes[indexSubPop].deathRate())
        {
            double carbonCReleased(0), carbonDReleased(0);
            deleteMicrobe(carbonCReleased,carbonDReleased,vaIndex);
            m_sizeC[positionM] += (1-m_erosion)*carbonCReleased;
            m_sizeD[positionM] += (1-m_erosion)*carbonDReleased;
            bool vide(true);
            for (int i1(0); i1<m_positions.size(); i1++)
            {
                for (int i2(0); i2<m_positions[i1].size(); i2++)
                    if (m_positions[i1][i2] == positionM) vide = false;
            }
            if (vide) m_types[positionM] = 0;
        }
    }


}


double Population::getphi(int const& indexSubPop) const
{
    return m_microbeTypes[indexSubPop].getphi();
}




void Population::calculusStepMeanVar(double const& timeStep, double &meanZ, double &varZ,
                                 double &meanC, double &varC, double &meanD,
                                 double &varD, double &meanM, double &varM,
                                 double &meanaccessDOC, double &varaccessDOC) const
{
    meanZ += timeStep*totalZ()/((double)m_length*m_length);
    meanC += timeStep*totalC()/((double)m_length*m_length);
    meanD += timeStep*totalD()/((double)m_length*m_length);
    meanM += timeStep*sizePopM()/((double)m_length*m_length);

    for (int k(0); k<m_length*m_length; k++)
    {
        varZ += timeStep*m_sizeZ[k]*m_sizeZ[k]/((double)m_length*m_length);
        varC += timeStep*m_sizeC[k]*m_sizeC[k]/((double)m_length*m_length);
        varD += timeStep*m_sizeD[k]*m_sizeD[k]/((double)m_length*m_length);
        varM += timeStep*m_sizeM[k]*m_sizeM[k]/((double)m_length*m_length);
        meanaccessDOC += timeStep*m_sizeD[k]*m_sizeM[k]/sizePopM();
        varaccessDOC += timeStep*m_sizeD[k]*m_sizeD[k]*m_sizeM[k]/sizePopM();
    }
}


void Population::calculusMeanVarBiomass(double &meanZ, double &varZ,
                                 double &meanC, double &varC, double &meanD,
                                 double &varD, double &meanM, double &varM,
                                 double &meanaccessDOC, double &varaccessDOC,
                                 double const& normalParamK, double const& weightEnzyme,
                                 double const& weightCarbonC, double const& weightCarbonD) const
{
    meanZ = totalZ()/((double)(m_length*m_length));
    meanC = totalC()/((double)(m_length*m_length));
    meanD = totalD()/((double)(m_length*m_length));
    meanM = sizePopM()/((double)(m_length*m_length));

    varZ = 0; 
    varC = 0; 
    varD = 0; 
    varM = 0; 
    meanaccessDOC = 0; 
    varaccessDOC = 0;

    for (int k(0); k<m_length*m_length; k++)
    {
        meanaccessDOC += (m_sizeD[k]*m_sizeM[k]/sizePopM());
    }

    for (int k(0); k<m_length*m_length; k++)
    {
        varZ += (m_sizeZ[k]-meanZ)*(m_sizeZ[k]-meanZ)/((double)(m_length*m_length));
        varC += (m_sizeC[k]-meanC)*(m_sizeC[k]-meanC)/((double)(m_length*m_length));
        varD += (m_sizeD[k]-meanD)*(m_sizeD[k]-meanD)/((double)(m_length*m_length));
        varM += (m_sizeM[k]-meanM)*(m_sizeM[k]-meanM)/((double)(m_length*m_length));
        varaccessDOC += (m_sizeD[k]-meanaccessDOC)*(m_sizeD[k]-meanaccessDOC)*m_sizeM[k]/sizePopM();
    }

    varZ *= normalParamK*weightEnzyme*normalParamK*weightEnzyme;
    varC *= normalParamK*weightCarbonC*normalParamK*weightCarbonC;
    varD *= normalParamK*weightCarbonD*normalParamK*weightCarbonD;
    varaccessDOC *= normalParamK*weightCarbonD*normalParamK*weightCarbonD;

    meanZ *= normalParamK*weightEnzyme;
    meanC *= normalParamK*weightCarbonC;
    meanD *= normalParamK*weightCarbonD;
    meanaccessDOC *= normalParamK*weightCarbonD;
}




void Population::writeFile(std::string nameFile, double const& timeIndex, double const& normalParamK,
                           double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD)
{
    ofstream flux(nameFile.c_str(),ios::app);
    flux << timeIndex << " ";
    for (int k(0); k<m_length*m_length; k++) flux << m_sizeZ[k]*normalParamK*weightEnzyme << " ";
    for (int k(0); k<m_length*m_length; k++) flux << m_sizeC[k]*normalParamK*weightCarbonC << " ";
    for (int k(0); k<m_length*m_length; k++) flux << m_sizeD[k]*normalParamK*weightCarbonD << " ";
    vector<int> sizeMTotal(m_length*m_length,0);
    for (int j(0); j<m_positions.size(); j++)
    {
        for (int i(0); i<m_positions[j].size(); i++) sizeMTotal[m_positions[j][i]] += 1;
    }
    for (int k(0); k<m_length*m_length; k++) flux << sizeMTotal[k] << " ";
    flux << endl;
}






void Population::writeFile2types(std::string nameFile, double const& timeIndex, double const& normalParamK,
                                 double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD,
                                 double const& phi1, double const& phi2)
{
    ofstream flux(nameFile.c_str(),ios::app);
    flux << timeIndex << " ";
    for (int k(0); k<m_length*m_length; k++) flux << m_sizeZ[k]*normalParamK*weightEnzyme << " ";
    for (int k(0); k<m_length*m_length; k++) flux << m_sizeC[k]*normalParamK*weightCarbonC << " ";
    for (int k(0); k<m_length*m_length; k++) flux << m_sizeD[k]*normalParamK*weightCarbonD << " ";
    vector<int> sizeMTotal1(m_length*m_length,0);
    vector<int> sizeMTotal2(m_length*m_length,0);
    if (nbrTypes() == 2)
    {
        for (int i(0); i<m_positions[0].size(); i++) sizeMTotal1[m_positions[0][i]] += 1;
        for (int i(0); i<m_positions[1].size(); i++) sizeMTotal2[m_positions[1][i]] += 1;
    }
    else if (nbrTypes() == 1)
    {
        if (getphi(0) == phi1) {for (int i(0); i<m_positions[0].size(); i++) sizeMTotal1[m_positions[0][i]] += 1;}
        if (getphi(0) == phi2) {for (int i(0); i<m_positions[0].size(); i++) sizeMTotal2[m_positions[0][i]] += 1;}
    }
    else if (nbrTypes() > 2)
    {
        cout << "nbrTypes() is neither 0, nor 1, nor 2: stop" << endl;
        abort();
    }

    for (int k(0); k<m_length*m_length; k++) flux << sizeMTotal1[k] << " ";
    for (int k(0); k<m_length*m_length; k++) flux << sizeMTotal2[k] << " ";
    flux << endl;

}







void Population::writeFile2typesTotalAccess(std::string nameFile, double const& timeIndex, double const& normalParamK,
                                            double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD,
                                            double const& phi1, double const& phi2)
{
    ofstream flux(nameFile.c_str(),ios::app);
    flux << timeIndex << " ";
    vector<int> sizeMTotal1(m_length*m_length,0);
    double popTotal1(0);
    vector<int> sizeMTotal2(m_length*m_length,0);
    double popTotal2(0);
    if (nbrTypes() == 2)
    {
        for (int i(0); i<m_positions[0].size(); i++) {sizeMTotal1[m_positions[0][i]] += 1; popTotal1 += 1;}
        for (int i(0); i<m_positions[1].size(); i++) {sizeMTotal2[m_positions[1][i]] += 1; popTotal2 += 1;}
    }
    else if (nbrTypes() == 1)
    {
        if (getphi(0) == phi1) {for (int i(0); i<m_positions[0].size(); i++) {sizeMTotal1[m_positions[0][i]] += 1; popTotal1 += 1;}}
        if (getphi(0) == phi2) {for (int i(0); i<m_positions[0].size(); i++) {sizeMTotal2[m_positions[0][i]] += 1; popTotal2 += 1;}}
    }
    else if (nbrTypes() > 2)
    {
        cout << "nbrTypes() is neither 0, nor 1, nor 2: stop" << endl;
        abort();
    }

    double sumZ(0), sumC(0), sumD(0), accessDOC1(0), accessDOC2(0), accessRelatif1(0), accessRelatif2(0);
    for (int k(0); k<m_length*m_length; k++)
    {
        sumZ += m_sizeZ[k];
        sumC += m_sizeC[k];
        sumD += m_sizeD[k];
        if (sizeMTotal1[k] >= 1)
        {
            accessDOC1 += m_sizeD[k]*sizeMTotal1[k];
            accessRelatif1 += (double)sizeMTotal1[k]/((double)(sizeMTotal1[k]+sizeMTotal2[k]))*m_sizeD[k];
        }
        if (sizeMTotal2[k] >= 1)
        {
            accessDOC2 += m_sizeD[k]*sizeMTotal2[k];
            accessRelatif2 += (double)sizeMTotal2[k]/((double)(sizeMTotal1[k]+sizeMTotal2[k]))*m_sizeD[k];
        }
    }
    if (popTotal1 > 0) {accessDOC1 /= popTotal1; accessRelatif1 /= popTotal1;}
    if (popTotal2 > 0) {accessDOC2 /= popTotal2; accessRelatif2 /= popTotal2;}

    flux << sumZ*normalParamK*weightEnzyme << " " << sumC*normalParamK*weightCarbonC << " ";
    flux << sumD*normalParamK*weightCarbonD << " ";
    flux << accessDOC1*normalParamK*weightCarbonD << " " << accessRelatif1*normalParamK*weightCarbonD << " ";
    flux << accessDOC2*normalParamK*weightCarbonD << " " << accessRelatif2*normalParamK*weightCarbonD << " ";
    flux << popTotal1 << " " << popTotal2;
    flux << endl;
}


void Population::writeTotalAccess(std::string nameFile, double const& timeIndex, double const& normalParamK,
                                  double const& weightEnzyme, double const& weightCarbonC, double const& weightCarbonD,
                                  double const& write1, double const& write2, double const& write3, double const& write4, 
                                  double const& write5, double const& write6, double const& write7, double const& write8, 
                                  double const& write9, double const& write10)
{
    ofstream flux(nameFile.c_str(),ios::app);
    flux << timeIndex << " ";

    if (nbrTypes() > 1)
    {
        cout << "nbrTypes() is neither 0, nor 1: stop" << endl;
        abort();
    }


    flux << totalZ()*normalParamK*weightEnzyme << " " << totalC()*normalParamK*weightCarbonC << " ";
    flux << totalD()*normalParamK*weightCarbonD << " " << sizePopM() << " " << write1 << " " << write2;
    flux << " " << write3 << " " << write4 << " " << write5 << " " << write6 << " " << write7;
    flux << " " << write8 << " " << write9 << " " << write10;
    flux << endl;
}

