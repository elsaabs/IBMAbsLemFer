#include <iostream>
#include "./code/Population.h"
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char* argv[])
{
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //parameters to specify

    cout << "argc = " << argc << endl;
    if (argc != 4)
        {
            cout << "3 arguments dans l'ordre: ID cluster min, ID cluster max, ID job max" << endl;
            return -1;
        }

    string nameFile("a");
    ifstream lecture; //pour la lecture
    string nameFileFinal("scoresCommuns.txt");
    ofstream fluxEcriture(nameFileFinal.c_str());

    for (int k(atoi(argv[1])); k <= atoi(argv[2]); k++)
    {
    for (int l(0); l <= atoi(argv[3]); l++)
    {

        ostringstream oss;
        oss << "./communs/scoresCommuns" << k << "_" << l << ".txt";
        nameFile = oss.str();
        lecture.open(nameFile.c_str(),ios::in);
        if(lecture)  // si l'ouverture a réussi
        {
            string contenu;  // déclaration d'une chaîne qui contiendra la ligne lue
            getline(lecture, contenu);  // on met dans "contenu" la ligne
            fluxEcriture << contenu;  // on affiche la ligne
            fluxEcriture << endl;
            lecture.close();  // on ferme le fichier
        }

    }
    }


    return 0;
}
