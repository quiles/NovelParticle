



#include "Snap.h"

#undef min
#undef max

#include "ParticleNet.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include <stdlib.h>
#include "Functions.h"


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
string fName, fNameCom="", newnetwork="";
string fileOfNames;
int maxSteps=10000;
bool dynamic=false, gen=false, simC=false;
int dim=3;
float epsilon;
float threshold=0.5;
bool fixedcom=false;
PParticleNet Model;


void SaveMeasures(const char *filename){
//    printf("nodes:%d  edges:%d\n", Graph->GetNodes(), Graph->GetEdges());
//
//    PUNGraph UGraph = TSnap::ConvertGraph<PParticleNet>(Model); // undirected version of the graph

    TIntFltH BtwH, EigH, PRankH, CcfH, CloseH, HubH, AuthH;
//    printf("Computing...\n");
//    printf("Treat graph as DIRECTED: ");
    //printf(" PageRank... ");
    TSnap::GetPageRank(Model, PRankH, 0.85);
    //printf(" Hubs&Authorities...");
    TSnap::GetHits(Model, HubH, AuthH);
//    printf("\nTreat graph as UNDIRECTED: ");
//    printf(" Eigenvector...");           TSnap::GetEigenVectorCentr(Model, EigH);
    //printf(" Clustering...");
    TSnap::GetNodeClustCf(Model, CcfH);
    //printf(" Betweenness (SLOW!)...");
    TSnap::GetBetweennessCentr(Model, BtwH, 1.0);
//    printf(" Constraint (SLOW!)...");    TNetConstraint<PParticleNet> NetC(Model, true);
//    printf(" Closeness (SLOW!)...");

    for (TParticleNet::TNodeI NI = Model->BegNI(); NI < Model->EndNI(); NI++) {
        const int NId = NI.GetId();
        CloseH.AddDat(NId, TSnap::GetClosenessCentr<PParticleNet>(Model, NId, false));
    }
    FILE *F = fopen(filename, "w");
//    fprintf(F,"#Network: %s\n", InFNm.CStr());
    fprintf(F,"%%Nodes: %d\tEdges: %d\n", Model->GetNodes(), Model->GetEdges());
    fprintf(F,"%%NodeId\tDegree\tCloseness\tBetweennes\tClusteringCoefficient\tPageRank\tHubScore\tAuthorityScore\n");
    for (TParticleNet::TNodeI NI = Model->BegNI(); NI < Model->EndNI(); NI++) {
        const int NId = NI.GetId();
        const double DegCentr = Model->GetNI(NId).GetDeg();
        const double CloCentr = CloseH.GetDat(NId);
        const double BtwCentr = BtwH.GetDat(NId);
//        const double EigCentr = EigH.GetDat(NId);
//        const double Constraint = NetC.GetNodeC(NId);
        const double ClustCf = CcfH.GetDat(NId);
        const double PgrCentr = PRankH.GetDat(NId);
        const double HubCentr = HubH.GetDat(NId);
        const double AuthCentr = AuthH.GetDat(NId);
        fprintf(F, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", NId,
                DegCentr, CloCentr, BtwCentr, ClustCf, PgrCentr, HubCentr, AuthCentr);
    }
    fclose(F);
}

float GetBetweenness(){
    TIntFltH BtwH;
    TSnap::GetBetweennessCentr(Model, BtwH, 1.0);
    float bet=0.0;
    for (TParticleNet::TNodeI NI = Model->BegNI(); NI < Model->EndNI(); NI++) {
        const int NId = NI.GetId();
        bet += BtwH.GetDat(NId);
    }
    return bet;
}


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
        cout << "\t\t-dim value ->  define the dimension of the particle system [default: 3]\n";
        cout << "\t\t-ltr value ->  Link threshold [default: 0.5]\n";
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



void ModelDynamic(int exec){
    string saveName;
    int it, st, count, itc, i;
    char out[256], ext[256];
    bool firstIt=true;
    int cInfo, maxInfo;
    float nmi;
    
    FILE *stream;
    sprintf(out,"output_k5.txt");
    stream = fopen(out, "w");
    fprintf(stream,"%% #file #transient #communities #max-com-size #infomap #centroid_error\n");

    steps = maxSteps;
    ifstream file (fileOfNames.c_str());
    count = 0;
    
    while  (file >> fName){
        if (firstIt){
            Model = new TParticleNet(dim);
            Model->LoadFromFileMat(fName.c_str());
            Model->SetModelParameters(alpha, beta, 1.0);
            sprintf(out,"time_0.par");
            saveName = fName;
            saveName.replace(fName.size()-3,3,out);
//            cout << "Initial state: " << saveName << endl;
            Model->SaveParticlePosition(saveName.c_str());
            firstIt = false;
        }
        else {
            Model->ReloadNetworkMat(fName.c_str());
        }
        cout << "Running model on: " << fName << endl;

        if (!fNameCom.empty()) Model->LoadComFile(fNameCom.c_str());

        Model->SetTag(count);
        it = Model->RunModelMat(steps,minDR,threshold, verbose);

//        for (i=0 ; i<100 ; i++){
//            Model->RunByStep();
//            sprintf(out,"%d.par",i);
//            saveName = fName;
//            saveName.replace(fName.size()-3,3,out);
//cout << saveName << endl;
//            Model->SaveParticlePosition(saveName.c_str());
//        }
//        Model->RunByStep();
//        sprintf(out,"%d.par",i);
//        saveName = fName;
//        saveName.replace(fName.size()-3,3,out);
//cout << saveName << endl;
//        itc = Model->CommunityDetection3();
//        Model->SaveParticlePosition(saveName.c_str());


        itc = 0;
        itc = Model->CommunityDetection3();
//        Model->CommunityDetectionDB(0.2);
        cout << "T: " << count << " ";
        cout << "Steps: " << it << " ";
//        cout << "StepsC: " << itc << " ";
//        cout << "#Com: " << Model->getNumCommunities() << " ";
//        cout << "CE: " << Model->printCentroidsError() << " ";
//        cout << "SizeC: " << Model->sizeLargeCom() << " ";
//        cout << "DiffX: " << Model->getDiffX() << " ";
//        cout << "Var: " << Model->getVar() << " ";
        cout << "FR: " << Model->getNormFR(0) << endl;

//        cInfo = Model->Infomap(maxInfo);
        cInfo = maxInfo = 0;

//        cout << endl << endl;

        if (!fNameCom.empty()) nmi = Model->NMI();
        else nmi = 0.0;

//        sprintf(ext,"t%d.par",count);
//        saveName = fName;
//        saveName.replace(fName.size()-3,3,ext);
//        Model->SaveParticlePosition(saveName.c_str());

//        saveName = fName;
//        saveName.replace(fName.size()-3,3,"for");
//        Model->SaveParticleForces(saveName.c_str());
//
//        saveName = fName;
//        saveName.replace(fName.size()-3,3,"mes");
//        SaveMeasures(saveName.c_str());        

        fprintf(stream,"%d %d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %d %d\n",
                count++, //1
                it, //2
                Model->getNetDegree(), //3
                Model->getNormFR(0), //4
                Model->getNormFR(1), //5
                Model->getNormFR(2), //6
                Model->getNormFR(3), //7
                Model->getVar(0), //8 (MSD distances from the centre of mass)
                Model->getVar(1), //9 (MSD distances from the centre of mass)
                Model->getVar(2), //10 (MSD distances from the centre of mass)
                Model->getVar(3), //11 (MSD distances from the centre of mass)
                Model->getDiffX(0), // 12 (RMSD of positions)
                Model->getDiffX(1), // 13 (RMSD of positions)
                Model->getDiffX(2), // 14 (RMSD of positions)
                Model->getDiffX(3), // 15 (RMSD of positions)
                Model->getDiffR(0), // 16 (RMSD of repulsion interactions)
                Model->getDiffR(1), // 17 (RMSD of repulsion interactions)
                Model->getDiffR(2), // 18 (RMSD of repulsion interactions)
                Model->getDiffR(3), // 19 (RMSD of repulsion interactions)
                Model->getNumCommunities(), // 20
                Model->sizeLargeCom()); // 21

        
        
    }
    fclose(stream);
}


/*
void ModelDynamic(){
    string saveName;
    int it, st, count, itc;
    char out[256];
    bool firstIt=true;
    int cInfo, maxInfo;
 
    FILE *stream;
    stream = fopen("output.txt", "w");
    fprintf(stream,"%% #file #transient #communities #max-com-size #infomap #centroid_error\n");

    steps = maxSteps;
    ifstream file (fileOfNames.c_str());
//cout << "[";        
count = 1;
    while  (file >> fName){
        if (firstIt){
            Model = new TParticleNet(dim);
            Model->LoadFromFile(fName.c_str());
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

        itc = Model->CommunityDetection3();
        cout << "Steps: " << it << " ";
        cout << "StepsC: " << itc << " ";
        cout << "#Com: " << Model->getNumCommunities() << " ";
        cout << "CE: " << Model->printCentroidsError() << " ";
        cout << "SizeC: " << Model->sizeLargeCom() << " ";
        cout << "FR: " << Model->getNormFR() << endl;

//        cInfo = Model->Infomap(maxInfo);
        cInfo = maxInfo = 0;

        cout << endl << endl;


        fprintf(stream,"%d %d %d %d %d %d %d %.5f\n",count++,it,itc,
                                            Model->getNumCommunities(),
                                            Model->sizeLargeCom(), 
                                            cInfo, 
                                            maxInfo,
                                            Model->printCentroidsError());
//cout << count++ << " " << it << " " << Model->getNumCommunities() << " " << Model->printCentroidsError() << "\n";
        saveName = fName;
        saveName.replace(fName.size()-3,3,"par");
        Model->SaveParticlePosition(saveName.c_str());

        saveName = fName;
        saveName.replace(fName.size()-3,3,"mes");
        SaveMeasures(saveName.c_str());

        
//        cout << "Current state file: " << saveName << endl;
        
//        saveName.replace(fName.size()-3,3,"com");
//        Model->SaveCommunities(saveName.c_str());
//        cout << "Current community structure: " << saveName << endl << endl;
    }
//cout << "]\n" << endl;
}
*/
void ModelByStep(){
    string saveName;
    char out[256];
    int i;
    clock_t ini,end;
    
    ini = clock();
    
    Model = new TParticleNet(dim);
    Model->LoadFromFile(fName.c_str());

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
        Model->RunModel(0, minDR, verbose);
//        Model->RunByStep();
        if (i>=0 && i%savestep==0) {
            sprintf(out,"time_%d.par",i);
            saveName = fName;
            saveName.replace(fName.size()-3,3,out);
            cout << "Step: " << i << " - " << saveName.c_str() << endl;
            Model->SaveParticlePosition(saveName.c_str());
        }
    }

    cout << "Detecting clusters...\n";
    Model->CommunityDetection3();
    end = clock();

    Model->CommunityDetectionDB(0.5);
    
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
    string saveName;
    char out[256];
    int i;
    clock_t ini,end;
    
    ini = clock();
    
    Model = new TParticleNet(dim);
//    Model->LoadFromFileMat(fName.c_str());
    Model->LoadFromFile(fName.c_str());
    
//    Model = TParticleNet::LoadFromFile(fName.c_str(),dim);
    if (!fNameCom.empty()) Model->LoadComFile(fNameCom.c_str());
    
    sprintf(out,"time_0.par");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveParticlePosition(saveName.c_str());
    cout << "Initial state: " << saveName << endl;
    
    Model->SetModelParameters(alpha, beta, 1.0);
    cout << "Model running...\n";

    Model->RunModel2(maxSteps, minDR, verbose);
//    Model->RunModelMat(maxSteps, minDR, threshold, verbose);

    cout << "Detecting clusters...\n";
    Model->CommunityDetection3();
    end = clock();
    
    sprintf(out,"cen");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveCentroids(saveName.c_str());
    
    cout << "# of communities detected: " << Model->getNumCommunities() << endl;
    cout << "FR norm: " << Model->getNormFR(0) << endl;

    if (!fNameCom.empty()) cout << "NMI: " << Model->NMI() << endl;
    cout << "Elapsed time (s): " << ((float)(end-ini))/CLOCKS_PER_SEC << endl;
    
    sprintf(out,"mes");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    SaveMeasures(saveName.c_str());

    if (!newnetwork.empty()) Model->SaveNetworkFromParticle(newnetwork.c_str(),epsilon);
    
    sprintf(out,"time_final.par");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveParticlePosition(saveName.c_str());
    cout << "Final state: " << saveName << endl;    

    FILE *st;
    st = fopen("Exp_dim.txt","a+");
    fprintf(st,"%d %d %.3f %.3f %.3f %.3f %.3f",
                Model->getNumCommunities(), //1
                Model->sizeLargeCom(),  //2
                Model->NMI(), //3
                Model->printCentroidsError(), //4
                Model->getNormFR(0), //5
                Model->getVar(0), //6
                Model->getDiffX(0)); //7

}

void ExpDim(float zout, int ite){
    string saveName;
    char out[256];
    int i;
    clock_t ini,end;
    int steps;
    
    ini = clock();
    
    Model = new TParticleNet(dim);
    Model->LoadFromFile(fName.c_str());
    
//    Model = TParticleNet::LoadFromFile(fName.c_str(),dim);
    if (!fNameCom.empty()) Model->LoadComFile(fNameCom.c_str());
  
    sprintf(out,"time_0.par");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveParticlePosition(saveName.c_str());
    cout << "Running ite: " << ite << endl;
    
    Model->SetModelParameters(alpha, beta, 1.0);
    steps = Model->RunModel(maxSteps, minDR, verbose);
    Model->CommunityDetection3();    

    FILE *st;
    sprintf(out,"output_%.2f_%d.txt",zout,dim);
    st = fopen(out,"a+");
    fprintf(st,"%.2f %d %d %d %.5f %.3f %.3f %.3f %.3f\n",
		zout,                
		steps,
		Model->getNumCommunities(), //2
                Model->sizeLargeCom(),  //3
                Model->NMI(), //4
                Model->printCentroidsError(), //5
                Model->getNormFR(0), //6
                Model->getVar(0), //7
                Model->getDiffX(0) //8
    );
    fclose(st);
}



//void ModelLote(){
//    PParticleNet Model;
//    char fcom[256], fnet[256];
//    int i;
//    clock_t ini,end;
//    float mu;
//
//    alpha = 1.0;
//    beta = 0.2;
//    steps = 100;
//    
//    for (mu=0.6 ; mu<=0.81 ; mu+=0.02){
//        ini = clock();
//
//        sprintf(fnet,"./data1000S/net_m%.2f_r1.dat",mu);
//        sprintf(fcom,"./data1000S/com_m%.2f_r1.dat",mu);
//        cout << "File: " << fnet << endl;
//
//        Model = TParticleNet::LoadFromFile(fnet,dim);
//        Model->LoadComFile(fcom);
//      
//        Model->SetModelParameters(alpha, beta, 1.0);
//        for (i=1 ; i<steps ; i++){
//            Model->RunByStep();
//        }
//        Model->CommunityDetection3();
//        end = clock();
//
//        cout << "# of communities detected: " << Model->getNumCommunities() << endl;
//        cout << "NMI: " << Model->NMI() << endl;
//        cout << "Elapsed time (s): " << ((float)(end-ini))/CLOCKS_PER_SEC << endl;
//    }
//}

void SimChange(){
    string saveName;
    char out[256];
    int i,count=0, it;
    FILE *stream;
    
    Model = new TParticleNet(dim);
    Model->LoadFromFile(fName.c_str());
    Model->SetModelParameters(alpha, beta, 1.0);

//    printf("Debug: %s\n",fName.c_str() );
//    return;
    
    stream = fopen("output.txt", "w");
    
    cout << "Model running...\n";
    
    for (i=1 ; i<100 ; i++){
        cout << i << endl;
        it = Model->RunModel(maxSteps, minDR, verbose);
        sprintf(out,"time_%d.par",i);
        saveName = fName;
        saveName.replace(fName.size()-3,3,out);
        Model->SaveParticlePosition(saveName.c_str());
        fprintf(stream,"%d %d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n",
                count++, //1
                it, //2
                Model->getNetDegree(), //3
                Model->getNormFR(0), //4
                Model->getNormFR(1), //5
                Model->getNormFR(2), //6
                Model->getNormFR(3), //7
                Model->getVar(0), //8 (MSD distances from the centre of mass)
                Model->getVar(1), //9 (MSD distances from the centre of mass)
                Model->getVar(2), //10 (MSD distances from the centre of mass)
                Model->getVar(3), //11 (MSD distances from the centre o f mass)
                Model->getDiffX(0), // 12 (RMSD of positions)
                Model->getDiffX(1), // 13 (RMSD of positions)
                Model->getDiffX(2), // 14 (RMSD of positions)
                Model->getDiffX(3), // 15 (RMSD of positions)
                Model->getDiffR(0), // 16 (RMSD of repulsion interactions)
                Model->getDiffR(1), // 17 (RMSD of repulsion interactions)
                Model->getDiffR(2), // 18 (RMSD of repulsion interactions)
                Model->getDiffR(3)); // 19 (RMSD of repulsion interactions)
        Model->ChangeNetwork(0,0.01);
    }
    fclose(stream);


    sprintf(out,"time_%d.par",i);
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveParticlePosition(saveName.c_str());
    cout << "Final state: " << saveName << endl;
   
}


int main(int argc,char *argv[]){
    int i=0;
    srand (time(NULL));
    bool byStep=false;
    float k_eps;
    bool mst=false;
    int opt_gen;
    float zout;
    
    
//    for (float r=-1.0 ; r<1.0 ; r+=0.01){
//        zout = r *M_PI/2.0;
//        printf("Tan(%.2f) = %.2f\n",r,1.0/tan(zout));
//    }
//    
//    return 0;
    
    printf("Running: \n$ ");
    for (i=0 ; i<argc ; i++) {
        printf("%s ",argv[i]);
    }
    printf("\n\n");
    
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
        else if (strcmp(argv[1],"-generate")==0){
            stream = fopen(argv[2], "r");
            if (!stream){ Message(1); return 0;}
            fclose(stream);
            fileOfNames = argv[2];
            gen = true;
        }
        else if (strcmp(argv[1],"-change")==0){
            stream = fopen(argv[2], "r");
            if (!stream){ Message(1); return 0;}
            fclose(stream);
            fName = argv[2];
            simC= true;
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
                fixedcom = true;
            }
            else if (strcmp(argv[i],"-max") == 0){
                if (++i>=argc) break;
                maxSteps = strtol(argv[i],NULL, 10);
            }
            else if (strcmp(argv[i],"-dim") == 0){
                if (++i>=argc) break;
                dim = strtol(argv[i],NULL, 10);
            }
            else if (strcmp(argv[i],"-tr") == 0){
                if (++i>=argc) break;
                minDR = (float)strtod(argv[i],NULL);
            }
            else if (strcmp(argv[i],"-ltr") == 0){
                if (++i>=argc) break;
                threshold = (float)strtod(argv[i],NULL);
            }
            else if (strcmp(argv[i],"-a") == 0){
                if (++i>=argc) break;
//                alpha = (float)strtod(argv[i],NULL);
		alpha = 1.0;
		zout = (float)strtod(argv[i],NULL);
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
            else if (strcmp(argv[i],"-eps") == 0){
                if (++i>=argc) break;
                k_eps = (float)strtod(argv[i],NULL);
                opt_gen = 1;
            }
            else if (strcmp(argv[i],"-knn") == 0){
                if (++i>=argc) break;
                k_eps = strtol(argv[i],NULL, 10);
                opt_gen = 2;
            }
            else if (strcmp(argv[i],"-mst") == 0){
                mst = true;
            }
        }

//        for (i=0 ; i<10 ; i++) {
//            ExpDim(zout,i);
//        }
//        return 0;
        
        if (simC) SimChange();
        else if (byStep) ModelByStep();
        else if (dynamic) ModelDynamic(0);
        else if (gen) GenerateNetworks(fileOfNames.c_str(), opt_gen, k_eps, mst);
        else {
            cout << "Model0\n";
            Model0();
        }
    }


    return 0;
}
