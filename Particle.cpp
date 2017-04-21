



#include "Snap.h"

#undef min
#undef max

#include "ParticleNet.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include <stdlib.h>

using namespace std;

void PrintGStats(const char s[], PUNGraph Graph) {
    cout << "Graph: " << s << " nodes: " << Graph->GetNodes()
    << " edges: " << Graph->GetEdges() << endl;
}

long int numCom=0;
float alpha=1.0;
float beta=0.30, betaM=1.0, betaS=0.1;
bool saveStates=false;
bool timeVarying=false;
bool snapFormat=false;
bool searchBeta=false;
bool verbose=false;
int steps=100;
int savestep=100;
float minDR=0.1;
string fName, fNameCom="";
string fileOfNames;
int maxSteps=1000;
bool dynamic=false;


void Message(int i){
    if (i==0){
        cout << "Community Detection - Particle's Model\n\n";
        cout << "\tUsage (1) : $ ./Particle filename.dat [options]\n\n";
        cout << "\tUsage (2) : $ ./Particle -dynamic filename.dat [options]\n\n";
        cout << "\t[options]: \n";
        cout << "\t\t-a value -> attractive parameter (alpha) [default: 1.0]\n";
        cout << "\t\t-b value -> repulsive parameter (beta) [default: 0.1]\n";
        cout << "\t\t-tr value ->  define the stop condition - theta_R  [default: 1.0]\n";
        cout << "\t\t-max value ->  define the maximum number of steps (iterations)  [default: 1000]\n";
        //        cout << "\t\t-ss number_of_steps -> save states (files.par & files.cen)\n";
        //        cout << "\t\t-sf -> load the input file using the snapFormat\n";
        cout << "\t\t-v -> verbose\n";
        cout << "\t\t-sb initial_beta final_beta step -> searching beta [available only for (1)]\n";
        cout << "\n\n";
        cout << "\tExample (1): $ ./Particle rede001.dat -b 0.3 -tr 1.0 - v\n\n";
        cout << "\tExample (2): $ ./Particle -dynamic fileOfNames.dat -b 0.3 \n";
        cout << "\t\tObs: fileOfNames.dat must contain the names of the network files, one file per line\n\n";
    }
    else if (i==1) cout << "error: no valid input files\n\n";
}



/*
 void ModelByStep(){
 PParticleNet Model;
 Model = TParticleNet::LoadFromFile("./net_sample.dat");
 string saveName;
 char out[256];
 int it, st, i;
 //    clock_t ini,end;
 //    ini = clock();
 fName = "./net_sample.dat";
 
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
 //        cout << "Step: " << i << " - " << saveName.c_str() << endl;
 Model->SaveParticlePosition(saveName.c_str());
 }
 }
 */

void ModelDynamic(){
    PParticleNet Model;
    string saveName;
    int it, st, count;
    char out[256];
    bool firstIt=true;
    
    FILE *stream;
    stream = fopen("output.txt", "w");
    
    ifstream file (fileOfNames.c_str());
//cout << "[";        
count = 1;
    while  (file >> fName){
        if (firstIt){
            Model = TParticleNet::LoadFromFile(fName.c_str());
            Model->SetModelParameters(alpha, beta, 1.0);
            sprintf(out,"time_0.par");
            saveName = fName;
            saveName.replace(fName.size()-3,3,out);
            cout << "Initial state: " << saveName << endl;
            Model->SaveParticlePosition(saveName.c_str());
            firstIt = false;
        }
        else {
            Model->ReloadNetwork(fName.c_str());
        }
        cout << "Running model on: " << fName << endl;

        it = Model->RunModel(steps,minDR,verbose);

        st = 0; //Model->CommunityDetection3();
        cout << "Total steps: " << it << endl;
//        cout << "# of communities detected: " << Model->getNumCommunities() << endl;
        cout << "Accumulated centroid error: " << Model->printCentroidsError() << endl;

        fprintf(stream,"%d %d %d %.5f\n",count++,it,Model->getNumCommunities(),Model->printCentroidsError());
//cout << count++ << " " << it << " " << Model->getNumCommunities() << " " << Model->printCentroidsError() << "\n";
        saveName = fName;
        saveName.replace(fName.size()-3,3,"par");
        Model->SaveParticlePosition(saveName.c_str());
//        cout << "Current state file: " << saveName << endl;
        
//        saveName.replace(fName.size()-3,3,"com");
//        Model->SaveCommunities(saveName.c_str());
//        cout << "Current community structure: " << saveName << endl << endl;
    }
//cout << "]\n" << endl;
}


void ModelByStep(){
    PParticleNet Model;
    string saveName;
    char out[256];
    int i;
    clock_t ini,end;
    
    ini = clock();
    
    Model = TParticleNet::LoadFromFile(fName.c_str());
    if (!fNameCom.empty()) Model->LoadComFile(fNameCom.c_str());
  
    sprintf(out,"time_0.par");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveParticlePosition(saveName.c_str());
    cout << "Initial state: " << saveName << endl;
    
    Model->SetModelParameters(alpha, beta, 1.0);
    cout << "Model running...\n";

//    Model->RunModel(steps, 0.1, false);

    for (i=1 ; i<steps ; i++){
        Model->RunByStep();
        if (i>=0 && i%savestep==0) {
            sprintf(out,"time_%d.par",i);
            saveName = fName;
            saveName.replace(fName.size()-3,3,out);
            cout << "Step: " << i << " - " << saveName.c_str() << endl;
            Model->SaveParticlePosition(saveName.c_str());
//            if (!fNameCom.empty()) {
//                Model->CommunityDetection3();
//                cout << "NMI: " << Model->NMI() << endl;
//            }
        }
    }

    cout << "Detecting clusters...\n";
    Model->CommunityDetection3();
    end = clock();
    
    cout << "# of communities detected: " << Model->getNumCommunities() << endl;
    if (!fNameCom.empty()) cout << "NMI: " << Model->NMI() << endl;
    cout << "Elapsed time (s): " << ((float)(end-ini))/CLOCKS_PER_SEC << endl;
    
    sprintf(out,"time_%d.par",i);
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveParticlePosition(saveName.c_str());
    cout << "Final state: " << saveName << endl;    
}


void Model0(){
    PParticleNet Model;
    string saveName;
    char out[256];
    int i;
    clock_t ini,end;
    
    ini = clock();
    
    Model = TParticleNet::LoadFromFile(fName.c_str());
    if (!fNameCom.empty()) Model->LoadComFile(fNameCom.c_str());
  
    sprintf(out,"time_0.par");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveParticlePosition(saveName.c_str());
    cout << "Initial state: " << saveName << endl;
    
    Model->SetModelParameters(alpha, beta, 1.0);
    cout << "Model running...\n";

    Model->RunModel(steps, minDR, verbose);

    cout << "Detecting clusters...\n";
    Model->CommunityDetection3();
    end = clock();
    
    cout << "# of communities detected: " << Model->getNumCommunities() << endl;
    if (!fNameCom.empty()) cout << "NMI: " << Model->NMI() << endl;
    cout << "Elapsed time (s): " << ((float)(end-ini))/CLOCKS_PER_SEC << endl;
    
    sprintf(out,"time_final.par");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveParticlePosition(saveName.c_str());
    cout << "Final state: " << saveName << endl;    
}


void ModelLote(){
    PParticleNet Model;
    char fcom[256], fnet[256];
    int i;
    clock_t ini,end;
    float mu;

    alpha = 1.0;
    beta = 0.2;
    steps = 100;
    
    for (mu=0.6 ; mu<=0.81 ; mu+=0.02){
        ini = clock();

        sprintf(fnet,"./data1000S/net_m%.2f_r1.dat",mu);
        sprintf(fcom,"./data1000S/com_m%.2f_r1.dat",mu);
        cout << "File: " << fnet << endl;

        Model = TParticleNet::LoadFromFile(fnet);
        Model->LoadComFile(fcom);
      
        Model->SetModelParameters(alpha, beta, 1.0);
        for (i=1 ; i<steps ; i++){
            Model->RunByStep();
        }
        Model->CommunityDetection3();
        end = clock();

        cout << "# of communities detected: " << Model->getNumCommunities() << endl;
        cout << "NMI: " << Model->NMI() << endl;
        cout << "Elapsed time (s): " << ((float)(end-ini))/CLOCKS_PER_SEC << endl;
    }
}


int main(int argc,char *argv[]){
    int i=0;
    srand (time(NULL));
    bool byStep=false;

//    ModelLote(); return 0;
    
    if (argc <= 1) Message(0);
    else {
        FILE *stream;
        if (strcmp(argv[1],"-dynamic")==0){
            stream = fopen(argv[2], "r");
            if (!stream){ Message(1); return 0;}
            fclose(stream);
            fileOfNames = argv[2];
            dynamic=true;
        }
        else {
            stream = fopen(argv[1], "r");
            if (!stream){ Message(1); return 0;}
            fclose(stream);
            fName = argv[1];
        }
        for (i=2 ; i<argc ; i++){
            if (strcmp(argv[i],"-c") == 0){
                if (++i>=argc) break;
                numCom = strtol(argv[i],NULL, 10);
            }
            else if (strcmp(argv[i],"-max") == 0){
                if (++i>=argc) break;
                maxSteps = strtol(argv[i],NULL, 10);
            }
            else if (strcmp(argv[i],"-tr") == 0){
                if (++i>=argc) break;
                minDR = (float)strtod(argv[i],NULL);
            }
            else if (strcmp(argv[i],"-a") == 0){
                if (++i>=argc) break;
                alpha = (float)strtod(argv[i],NULL);
            }
            else if (strcmp(argv[i],"-cf") == 0){
                if (++i>=argc) break;
                fNameCom = argv[i];
            }
            else if (strcmp(argv[i],"-b") == 0){
                if (++i>=argc) break;
                beta = (float) strtod(argv[i],NULL);
            }
            else if (strcmp(argv[i],"-ss") == 0){
                saveStates = true;
                if (++i>=argc) break;
                savestep = strtol(argv[i],NULL,10);
                if (savestep == 0) --i;
            }
            else if (strcmp(argv[i],"-ms") == 0){
                if (++i>=argc) break;
                steps = strtol(argv[i],NULL,10);
                if (steps == 0) --i;
            }
            else if (strcmp(argv[i],"-tv") == 0){
                timeVarying = true;
            }
            else if (strcmp(argv[i],"-sf") == 0){
                snapFormat = true;
            }
            else if (strcmp(argv[i],"-v") == 0){
                verbose = true;
            }
            else if (strcmp(argv[i],"-bs") == 0){
                byStep = true;
            }
            //            else if (strcmp(argv[i],"-dynamic")==0){
            //                if (++i>=argc) break;
            //                fileOfNames = argv[i];
            //                dynamic = true;
            //                cout << "Dynamic: " << fileOfNames << endl;
            //            }
            else if (strcmp(argv[i],"-sb") == 0){
                searchBeta = true;
                if (++i>=argc) break;
                beta = (float) strtod(argv[i],NULL);
                if (++i>=argc) break;
                betaM = (float) strtod(argv[i],NULL);
                if (++i>=argc) break;
                betaS = (float) strtod(argv[i],NULL);
            }
        }

        if (byStep) ModelByStep();
//        else if (searchBeta) ModelSearchBeta();
        else if (dynamic) ModelDynamic();
        else Model0();
    }


    return 0;
}
