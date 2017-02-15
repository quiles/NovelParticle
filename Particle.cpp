
#include "Snap.h"

#undef min
#undef max

#include "ParticleNet.h"
#include <iostream>
//#include <fstream>
//#include <time.h>
//#include <string>
//#include <stdlib.h>

using namespace std;

long int numCom=0;
float alpha=1.0;
float beta=0.50, betaM=1.0, betaS=0.1;
bool saveStates=false;
bool timeVarying=false;
bool snapFormat=false;
bool searchBeta=false;
bool verbose=false;
int steps=1000;
float minDR=1.0;
string fName, fNameCom;
string fileOfNames;
int maxSteps=1000;
bool dynamic=false;



void ModelByStep(){
    PParticleNet Model;
    Model = TParticleNet::LoadFromFile("./netGN_sample.dat");
    string saveName;
    char out[256];
    int it, st, i;
//    clock_t ini,end;    
//    ini = clock();
    fName = "./netGN_sample.dat";
    
    sprintf(out,"time_0.par");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveParticlePosition(saveName.c_str());
    cout << "Initial state: " << saveName << endl;
    
//    Model->SetModelParameters(alpha, beta, 1.0);
    cout << "Model running...\n";
    for (i=0 ; i<100 ; i++){
        Model->RunByStep();
        sprintf(out,"time_%d.par",i);
        saveName = fName;
        saveName.replace(fName.size()-3,3,out);
        cout << "Step: " << i << " - " << saveName.c_str() << endl;
        Model->SaveParticlePosition(saveName.c_str());
    }
}


int main(int argc,char *argv[]){
    ModelByStep();

    return 0;
}
