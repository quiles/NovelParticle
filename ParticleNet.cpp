

#include "ParticleNet.h"

#include <time.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <limits>

#define DT 0.1
#define AGE 2


TNode::~ TNode(){

}

TParticleNet::TParticleNet(int dim){
    TCentroid centroid;
    PDIM = dim;
    for (i=0 ; i<PDIM ; i++) centroid.x[i] = 0;
    centroid.comm_id = 1;
    addCentroidThreshold = 0.5;
    numCommunities = 0;
    nextComId = 2; // id of the next detected community (used to identify the centroids/communities);
    Centroids.push_back(centroid);
    oldRR2 = 0.0;
    oldRR = 0.0;
    RR = 0.0;
    numClusters = 0;
}

void TParticleNet::LoadComFile(const char *filename){
    FILE *stream;
    PParticle particle;
    TParticleNet::TNodeI NI;
    int id, com;
    int maxCom=0;

    stream = fopen(filename, "r+");
    if (!stream) return;
    
    while (fscanf(stream, "%d %d", &id, &com) == 2){
        com++;
        if (IsNode(id)) {
            particle = GetNDat(id);
            particle->indexReal = com;
            if (com>maxCom) maxCom = com;
        }
        
//        cout << "ID: " << id << " - " << com << endl;
    }
    numCommunities = maxCom;

}


void TParticleNet::LoadFromFile(const char *filename){
    FILE *stream;
    stream = fopen(filename, "r+");
    if (!stream) {
        cout << "ERROR: File not found!\n";
        return;
    }
    int i;

    PParticle particle;
    PLinkData linkD;
    int id1, id2;

    while (fscanf(stream, "%d %d", &id1, &id2) == 2){
        if (!IsNode(id1)) {
            particle = new TParticle(PDIM);
            AddNode(id1,particle);
        }
        if (!IsNode(id2)) {
            particle = new TParticle(PDIM);
            AddNode(id2,particle);
        }
        if (!IsEdge(id1,id2)) {
            linkD = new TLinkData();
            linkD->age = AGE;
            AddEdge(id1, id2, linkD);
            AddEdge(id2, id1, linkD);
//            if (id1<id2) AddEdge(id1, id2, linkD);
//            else AddEdge(id2, id1, linkD);
        }
    }
    fclose(stream);
}

void TParticleNet::LoadFromFile2(const char *filename){
    FILE *stream;
    char tag[10];
    stream = fopen(filename, "r+");
    if (!stream) {
        cout << "ERROR: File not found!\n";
        return;
    }
    int i;
    
    PParticle particle;
    PLinkData linkD;
    int id1, id2;
    float weight;
    int parsing;
    bool nodes=true;
 
//    cout << "Nodes: ";
    while (fscanf(stream, "%s ",&tag) == 1){
        cout << tag << endl;
        if (tag[0] == 'N') {
            fscanf(stream, "%d", &id1);
            if (!IsNode(id1)) {
                particle = new TParticle(PDIM);
                AddNode(id1,particle);
            }
            cout << id1 << endl;
//            cout << id1 << " ";
        }
        else{
            if (nodes) cout << "\n\nLinks: \n";
            nodes = false;
            linkD = new TLinkData();
            parsing = fscanf(stream, "%d %d %f", &id1, &id2, &weight);
            if (parsing == 3) linkD->weight = weight;
            cout << id1 << " - " << id2 << " : " << linkD->weight << endl;
            linkD->age = AGE;
            AddEdge(id1, id2, linkD);
            AddEdge(id2, id1, linkD);
//            cout << "\t" << id1 << " --> " << id2 << endl;
        }
    }
    
    fclose(stream);
}

void TParticleNet::LoadFromFileMat(const char *filename){
    FILE *stream;
    char tag[10];
    stream = fopen(filename, "r+");
    if (!stream) {
        cout << "ERROR: File not found!\n";
        return;
    }
    int i;
    
    PParticle particle;
    PLinkData linkD;
    int id1, id2;
    float weight;
    int N;
    
    fscanf(stream,"%d", &N);
    for (id1=1; id1<=N ; id1++){
        particle = new TParticle(PDIM);
        AddNode(id1,particle);
    }
    for (id1=1 ; id1<N ; id1++){
        for (id2=id1+1 ; id2<=N ; id2++){
            fscanf(stream,"%f", &weight);
            if (weight>0.01){
                linkD = new TLinkData();
                linkD->weight = weight;
                linkD->age = AGE;
                AddEdge(id1, id2, linkD);
                AddEdge(id2, id1, linkD);
            }
        }
    }
    
    fclose(stream);
}


void TParticleNet::ChangeNetwork(int op, float per){
    TParticleNet::TNodeI NI;
    TParticleNet::TEdgeI EI;
    PLinkData linkD;
    int id1, id2, N;
    float guess;
    vector<int> listNodes;

    N = GetNodes();
    
    for (NI=BegNI(); NI<EndNI(); NI++) listNodes.push_back(NI.GetId());

cout << "\n\nChanges\n";
    
    switch(op){
        case 0: { // changes links randomly
            for (EI=BegEI() ; EI < EndEI() ; EI++){
                guess = (float)(rand() % 1001) / 1000.0;
//                cout << "Guess: " << guess << " -> ";
                id1 = EI.GetSrcNId();
                id2 = EI.GetDstNId();
//                cout << "Old Link: [" << id1 << ":" << id2 << "] ";
                if (guess <= per) {
                    // remove link
                    id1 = EI.GetSrcNId();
                    id2 = EI.GetDstNId();
//                    cout << "Old Link: [" << id1 << ":" << id2 << "] -> ";
                    linkD = GetEDat(id1,id2);
                    delete linkD;
                    DelEdge(id1,id2,false);
                    DelEdge(id2,id1,false);
                    // add new link
                    do {
                        id1 = listNodes[rand()%N];
                        id2 = listNodes[rand()%N];
                    } while (IsEdge(id1,id2) || id1==id2);
//                    cout << " -> New Link: [" << id1 << ":" << id2 << "]";
                    linkD = new TLinkData();
                    linkD->refreshed = true;
                    linkD->age = AGE;
                    AddEdge(id1, id2, linkD);
                    AddEdge(id2, id1, linkD);
                }
//                cout << endl;
            }
        }; break;
        case 1: { // add links  
            
        }; break;
        case 3: { // remove links
            
        }; break;
        case 4: { //
        }; break;
    }
//    cout << endl << endl;
//    for (EI=BegEI() ; EI < EndEI() ; EI++){
//        id1 = EI.GetSrcNId();
//        id2 = EI.GetDstNId();
//        cout << "Final Links: [" << id1 << ":" << id2 << "] \n";
//    }
//    cout << endl;
    
}

void TParticleNet::ReloadNetwork(const char *filename){
    FILE *stream;
    stream = fopen(filename, "r+");
    if (!stream) return;

    int id1, id2;
    PParticle particle;


// Inserir dados das arestas ao atualizar o grafo
// Verificar a condicao de remocao dos dados das arestas e nos da rede

    PLinkData linkD;
    TParticleNet::TNodeI NI;
    TParticleNet::TEdgeI EI;
    int idt;

    
    for (NI=BegNI(); NI<EndNI(); NI++) {
        particle = GetNDat(NI.GetId());
        particle->refreshed = false;
    }
    
    for (EI=BegEI() ; EI < EndEI() ; EI++){
        id1 = EI.GetSrcNId();
        id2 = EI.GetDstNId();
        linkD = GetEDat(id1,id2);
        linkD->refreshed = false;
    }
    

    while (fscanf(stream, "%d %d", &id1, &id2) == 2){
        if (!IsNode(id1)) {
            particle = new TParticle(PDIM);
            AddNode(id1,particle);
        }
        else {
            particle = GetNDat(id1);
            particle->refreshed = true;
        }

        if (!IsNode(id2)) {
            particle = new TParticle(PDIM);
            AddNode(id2,particle);
        }
        else {
            particle = GetNDat(id2);
            particle->refreshed = true;
        }

        if (!IsEdge(id1,id2)) {
            linkD = new TLinkData();
            linkD->refreshed = true;
            linkD->age = AGE;
            AddEdge(id1, id2, linkD);
            AddEdge(id2, id1, linkD);  
        }
        else {
            linkD = GetEDat(id1,id2);
            linkD->refreshed = true;
            linkD->age = AGE;
        }
    }
    fclose(stream);

    // remove unrefreshed links
    for (EI=BegEI() ; EI < EndEI() ; EI++){
        id1 = EI.GetSrcNId();
        id2 = EI.GetDstNId();
        linkD = GetEDat(id1,id2);
        if (linkD->refreshed) continue;
	    delete linkD;
        DelEdge(id1,id2,false);
        DelEdge(id2,id1,false);
    }


    // remove unrefreshed nodes
    for (NI=BegNI(); NI<EndNI(); NI++) {
//        cout << NI.GetDeg() << " ";
        particle = GetNDat(NI.GetId());
	if (!particle->refreshed) {
            delete particle;
            DelNode(NI.GetId());
        }
    }
//    cout << endl;
}

void TParticleNet::ReloadNetworkMat(const char *filename){
    FILE *stream;
    stream = fopen(filename, "r+");
    if (!stream) return;
    
    int id1, id2;
    PParticle particle;
    
    
    // Inserir dados das arestas ao atualizar o grafo
    // Verificar a condicao de remocao dos dados das arestas e nos da rede
    
    PLinkData linkD;
    TParticleNet::TNodeI NI;
    TParticleNet::TEdgeI EI;
    int idt, N;
    float weight;
    
    fscanf(stream, "%d", &N);
    for (id1=1 ; id1<N ; id1++){
        for (id2=id1+1 ; id2<=N ; id2++){
            fscanf(stream,"%f", &weight);
            if (weight){
                if (!IsEdge(id1,id2)) {
                    linkD = new TLinkData();
                    linkD->weight = weight;
                    linkD->age = AGE;
                    AddEdge(id1, id2, linkD);
                    AddEdge(id2, id1, linkD);
                }
                else {
                    linkD = GetEDat(id1,id2);
                    linkD->weight = weight;
                }
            }
            else {
                if (IsEdge(id1,id2)) {
                    linkD = GetEDat(id1,id2);
                    delete linkD;
                    DelEdge(id1,id2,false);
                    DelEdge(id2,id1,false);
                }
            }
        }
    }

    fclose(stream);
}


void TParticleNet::NewNode(int node_id){
/*
    PParticle particle;
    if (IsNode(node_id)) return;

    particle = new TParticle(PDIM);
    for (i=0 ; i<PDIM ; i++) particle->x[i] = (float)(rand()%2000 - 1000) / 10000.0;
    particle->index = NULL;
    particle->indexReal = 0;
    particle->cluster_id = 0;
    AddNode(node_id,particle);
*/
}

void TParticleNet::NewLink(int i, int j){
/*
    if (IsNode(i) && IsNode(j))
       if (!IsEdge(i,j)) {
           AddEdge(i,j);
//           AddEdge(j,i);
       }
*/
}

void TParticleNet::DeleteNode(int node_id){
/*
    PParticle particle;
    if (IsNode(node_id)) {
        particle = GetNDat(node_id);
        delete particle;
        DelNode(node_id);
    }
*/
}

void TParticleNet::DeleteLink(int i, int j){
/*
    PLinkData link;
    if (IsNode(i) && IsNode(j)){
       if (IsEdge(i,j)) {
           link = GetEDat(i,j);
           delete link;
           DelEdge(i,j,false);
           DelEdge(j,i,false);
       }
    }
*/
}



void TParticleNet::ResetParticles() {
    TCentroid centroid;
    PParticle particle;
    TParticleNet::TNodeI NI, NN;

    Centroids.clear();
    
    for (NI=BegNI(); NI<EndNI(); NI++) {
        particle = GetNDat(NI.GetId());
        for (i=0 ; i<PDIM ; i++) particle->x[i] = (float)(rand()%2000 - 1000) / 10000.0;
        particle->index = NULL;
        particle->indexReal = 0;
    }
    
    for (i=0 ; i<PDIM ; i++) centroid.x[i] = 0.0;
    centroid.comm_id = 1;
    nextComId = 2; // id of the next detected community (used to identify the centroids/communities);
    Centroids.push_back(centroid);    
}


//void TParticleNet::RunByStep(){
//    TParticleNet::TNodeI NI, NN;
//    TParticleNet::TEdgeI EI;
//    PParticle data1, data2;
//    PLinkData linkD;
//    
//    for (NI=BegNI(); NI<EndNI(); NI++) {
//        data1 = GetNDat(NI.GetId());
//        data1->ComCentrality = 0.0;
//        for (i=0 ; i<PDIM ; i++) {
//            data1->dA[i] = 0;
//            data1->dR[i] = 0;
//        }
//    }
//    
//    
//    // Attraction
//    for (EI=BegEI() ; EI < EndEI() ; EI++){
//        id1 = EI.GetSrcNId();
//        id2 = EI.GetDstNId();
//        // Loading node's data;
//        data1 = GetNDat(id1);
//        data2 = GetNDat(id2);
//        for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - data2->x[i];
//        r = 0;
//        for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
//
//        data1->ComCentrality += r; //pow(r,2.0);
//        data2->ComCentrality += r; //pow(r,2.0);
//
//        r = sqrt(r);
//
//        linkD = GetEDat(id1,id2);
//        linkD->distance = r;
//
//        if (r<1.0) r=1.0;
//
//        for (i=0 ; i<PDIM ; i++) sum[i] = diff[i]/r;
//        for (i=0 ; i<PDIM ; i++) {
//            data1->dA[i] -= sum[i];
//            data2->dA[i] += sum[i];
//        }
//    }
//
//    
//    // repulsion
//    for (NI=BegNI() ; NI < EndNI() ; NI++){
//        data1 = GetNDat(NI.GetId());
//        for (NN=NI ; NN < EndNI() ; NN++){
//            if (NI.IsNbrNId(NN.GetId()) || NI==NN) continue;
//            data2 = GetNDat(NN.GetId());         
//            for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - data2->x[i];
//            r = 0;
//            for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
//            r = sqrt(r);
//if (r<0.01) r = 0.01;
//            for (i=0 ; i<PDIM ; i++) sum[i] = exp(-r)*diff[i]/(r);
//
//            for (i=0 ; i<PDIM ; i++) {
//                data1->dR[i] -= sum[i];
//                data2->dR[i] += sum[i];
//            }
//        }
//    }
//        
//    for (NI=BegNI(); NI<EndNI(); NI++) {
//        if ((degree = NI.GetDeg())){
//            data1 = GetNDat(NI.GetId());            
//            for (i=0 ; i<PDIM ; i++) {
////                data1->dA[i] = (alpha*data1->dA[i])/(float)degree;
////                data1->dR[i] = (beta*data1->dR[i])/(float)degree;
////                data1->x[i] += DT*((data1->dA[i] - data1->dR[i]));
//                data1->x[i] += ((alpha*data1->dA[i] - beta*data1->dR[i]))/degree;
//            }
//        }
//    }
//}


void TParticleNet::SetTag(int t){
    tag = t;
}

int TParticleNet::RunModel(int maxIT, float minDR, bool verbose){
    TParticleNet::TNodeI NI, NN;
    TParticleNet::TEdgeI EI;
    PParticle data1, data2;
    int steps=0;
    float value;
    float dg1, dg2;
    float tempR;
    float diffX[4]={0.0,0.0,0.0,0.0}, dTemp[4]={0.0,0.0,0.0,0.0}, diffR[4]={0.0,0.0,0.0,0.0}, dTempR[4]={0.0,0.0,0.0,0.0};
    float Massa[PDIM];
    float N, Media;
    PLinkData linkD;

    N = GetNodes();
    RR = 0.0;

//    if (verbose) cout << " Step#   \\theta_R   \\Delta\\theta_R \n";

    netDegree = 0.0;
    for (NI=BegNI(); NI<EndNI(); NI++) {
        data1 = GetNDat(NI.GetId());
        data1->degree = NI.GetOutDeg();
        netDegree += data1->degree;
//        data1->meanDistance = (float)(data1->degree);
        for (i=0 ; i<PDIM ; i++) {
            data1->xO[i] = data1->x[i];
            data1->dRO[i] = data1->dR[i];
        }
    }
    netDegree /= N;

    
    do {
        
        for (NI=BegNI(); NI<EndNI(); NI++) {
            data1 = GetNDat(NI.GetId());
            data1->DistF = 0;
            for (i=0 ; i<PDIM ; i++) {
                data1->dA[i] = 0;
                data1->dR[i] = 0;
                Massa[i] = 0.0;
//                data1->R = data1->meanDistance / (float)(data1->degree);
//                data1->meanDistance = 0.0;
            }
        }

// Algorithm CORE (O(n^2))
        // Attraction O(m))
//        float meanD1, meanD2;
        float signal;
        for (EI=BegEI() ; EI < EndEI() ; EI++){
            id1 = EI.GetSrcNId();
            id2 = EI.GetDstNId();
            // Loading node's data;
            data1 = GetNDat(id1);
            data2 = GetNDat(id2);
            
            for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - data2->x[i];
            r = 0;
            for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
            r = sqrt(r);
            
            linkD = GetEDat(id1,id2);
            linkD->distance = r;

//            if (r<1) r=1.0;
            
            float teste;
            
            teste = exp(-0.5*r)+1.0;
//            teste = 1.0;
            
            for (i=0 ; i<PDIM ; i++) sum[i] = diff[i]/r;
            for (i=0 ; i<PDIM ; i++) {
                data1->dA[i] -= teste*sum[i];
                data2->dA[i] += teste*sum[i];
//                data1->dA[i] -= sum[i]*dg2;
//                data2->dA[i] += sum[i]*dg1;
            }
        }

        // repulsion O(n^2) ... n-k \approx n p/ k small
        for (NI=BegNI() ; NI < EndNI() ; NI ++){
            data1 = GetNDat(NI.GetId());
            if (data1->degree < 2) continue; // avoid repulsion on nodes with degree 1 or zero
            for (NN=NI ; NN < EndNI() ; NN++){
                if (NI.IsNbrNId(NN.GetId()) || NI==NN) continue;
                data2 = GetNDat(NN.GetId());
                if (data2->degree < 2) continue; // avoid repulsion on nodes with degree 1 or zero
                for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - data2->x[i];
                r = 0;
                for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
                r = sqrt(r);
                data1->DistF += r;
                data2->DistF += r;
                if (r<0.001) r = 0.001; // to avoid division per zero when the particles are in the same position
                for (i=0 ; i<PDIM ; i++) sum[i] = exp(-r)*diff[i]/(r);


//                float d1, d2;
                for (i=0 ; i<PDIM ; i++) {
//                    d1 = data1->degree;
//                    d2 = data2->degree;
//                    data1->dR[i] -= sum[i]*d2*d1;
//                    data2->dR[i] += sum[i]*d1*d2;
                    data1->dR[i] -= sum[i];
                    data2->dR[i] += sum[i];
                }
            }
        }

        ++steps;
        
        
        FILE * st;
        char str[256];
        if (steps%10 == 1) {
            sprintf(str,"Simulation.time.%d",steps);
            SaveParticlePosition(str);
        }


        oldRR = RR;
        RR = 0.0;
        for (NI=BegNI(); NI<EndNI(); NI++) {
            if ((degree = NI.GetOutDeg())){
                data1 = GetNDat(NI.GetId());
                tempR = 0;
                for (i=0 ; i<PDIM ; i++) {
//                    data1->xOL[i] = data1->x[i];
                    value = DT*((alpha*data1->dA[i] - beta*data1->dR[i]))/degree;
                    data1->x[i] += value;
                    data1->dAO[i] = value;
                    tempR += (data1->dR[i])*(data1->dR[i]);
                }
                RR += sqrt(tempR); //degree;
            }
        }
// end Algorithm CORE

        value = fabs(oldRR-RR); // / (float)Particles.size();
        if (verbose) cout << "Running: " << steps << " --> " << value << "  Measure: "  << endl; //" DiffX " << diffX << endl;
//        if (value < minDR) {
//            if (steps%2) maxIT = steps;
//            else maxIT = steps+1;
//        }
        
    } while (steps<maxIT);
    
    
    
    return steps;
}

int TParticleNet::RunModel2(int maxIT, float minDR, bool verbose){
    TParticleNet::TNodeI NI, NN;
    TParticleNet::TEdgeI EI;
    PParticle data1, data2;
    int steps=0;
    float value;
    float dg1, dg2;
    float tempR;
    float N, Media;
    PLinkData linkD;
    
    
    FILE *energy;
    
    energy = fopen("Energy.txt", "w");
    
    N = GetNodes();
    RR = 0.0;
    
    //    if (verbose) cout << " Step#   \\theta_R   \\Delta\\theta_R \n";
    
    netDegree = 0.0;
    for (NI=BegNI(); NI<EndNI(); NI++) {
        data1 = GetNDat(NI.GetId());
        data1->degree = NI.GetOutDeg();
        netDegree += data1->degree;
        //        data1->meanDistance = (float)(data1->degree);
        for (i=0 ; i<PDIM ; i++) {
            data1->xO[i] = data1->x[i];
            data1->dRO[i] = data1->dR[i];
            data1->vel[i] = 0.0;
        }
    }
    netDegree /= N;
    
    
    do {
        
        for (NI=BegNI(); NI<EndNI(); NI++) {
            data1 = GetNDat(NI.GetId());
            data1->DistF = 0;
            for (i=0 ; i<PDIM ; i++) {
                data1->dA[i] = 0;
                data1->dR[i] = 0;
//                Massa[i] = 0.0;
                //                data1->R = data1->meanDistance / (float)(data1->degree);
                //                data1->meanDistance = 0.0;
            }
        }
        
        // Algorithm CORE (O(n^2))
        // Attraction O(m))
        //        float meanD1, meanD2;
        float signal;
        for (EI=BegEI() ; EI < EndEI() ; EI++){
            id1 = EI.GetSrcNId();
            id2 = EI.GetDstNId();
            // Loading node's data;
            data1 = GetNDat(id1);
            data2 = GetNDat(id2);
            
            for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - data2->x[i];
            r = 0;
            for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
            r = sqrt(r);
            
            linkD = GetEDat(id1,id2);
            linkD->distance = r;
            
            //            if (r<1) r=1.0;
            
            float teste;
            
            teste = exp(-0.5*r)+1.0;
            //            teste = 1.0;
            
            for (i=0 ; i<PDIM ; i++) sum[i] = diff[i]/r;
            for (i=0 ; i<PDIM ; i++) {
                data1->dA[i] -= teste*sum[i];
                data2->dA[i] += teste*sum[i];
                //                data1->dA[i] -= sum[i]*dg2;
                //                data2->dA[i] += sum[i]*dg1;
            }
        }
        
        // repulsion O(n^2) ... n-k \approx n p/ k small
        for (NI=BegNI() ; NI < EndNI() ; NI ++){
            data1 = GetNDat(NI.GetId());
            if (data1->degree < 2) continue; // avoid repulsion on nodes with degree 1 or zero
            for (NN=NI ; NN < EndNI() ; NN++){
                if (NI.IsNbrNId(NN.GetId()) || NI==NN) continue;
                data2 = GetNDat(NN.GetId());
                if (data2->degree < 2) continue; // avoid repulsion on nodes with degree 1 or zero
                for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - data2->x[i];
                r = 0;
                for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
                r = sqrt(r);
                data1->DistF += r;
                data2->DistF += r;
                if (r<0.001) r = 0.001; // to avoid division per zero when the particles are in the same position
                for (i=0 ; i<PDIM ; i++) sum[i] = exp(-r)*diff[i]/(r);
                
                
                //                float d1, d2;
                for (i=0 ; i<PDIM ; i++) {
                    //                    d1 = data1->degree;
                    //                    d2 = data2->degree;
                    //                    data1->dR[i] -= sum[i]*d2*d1;
                    //                    data2->dR[i] += sum[i]*d1*d2;
                    data1->dR[i] -= sum[i];
                    data2->dR[i] += sum[i];
                }
            }
        }
        
        ++steps;
        
        
        FILE * st;
        char str[256];
        if (steps%10 == 1) {
            sprintf(str,"Simulation.time.%d",steps);
            SaveParticlePosition(str);
        }
        
        
        oldRR = RR;
        RR = 0.0;
        float CE=0.0, PE=0.0;
        float Vel=0.0, Force=0.0;
        
        for (NI=BegNI(); NI<EndNI(); NI++) {
            if ((degree = NI.GetOutDeg())){
                data1 = GetNDat(NI.GetId());
                tempR = 0;
                Vel = 0.0;
                Force = 0.0;
                for (i=0 ; i<PDIM ; i++) {
                    //                    data1->xOL[i] = data1->x[i];
                    value = ((alpha*data1->dA[i] - beta*data1->dR[i]) - 1.0*data1->vel[i]);
                    
                    Force += value*value;
                    
                    value /= (float)degree;
                    data1->vel[i] += value*DT;

                    Vel += data1->vel[i]*data1->vel[i];
                    
                    data1->x[i] += data1->vel[i]*DT;
                    data1->dAO[i] = value;
                    tempR += (data1->dR[i])*(data1->dR[i]);
                }
                CE += degree*sqrt(Vel);
                PE += sqrt(Force);
                
                RR += sqrt(tempR); //degree;
            }
        }
        // end Algorithm CORE
        
        fprintf(energy, "%d %.5f %.5f\n", steps, CE, PE);
        value = fabs(oldRR-RR); // / (float)Particles.size();
        if (verbose) cout << "Running: " << steps << " --> " << value << "  CE: " << CE << " PE: " << PE << endl; //" DiffX " << diffX << endl;
        //        if (value < minDR) {
        //            if (steps%2) maxIT = steps;
        //            else maxIT = steps+1;
        //        }
        
    } while (steps<maxIT);
    
    fclose(energy);
    
    return steps;
}


int TParticleNet::RunModelMat(int maxIT, float minDR, float threshold, bool verbose){
    TParticleNet::TNodeI NI, NN;
    TParticleNet::TEdgeI EI;
    PParticle data1, data2;
    int id1, id2, d1, d2;
    int steps=0;
    float value;
    float dg1, dg2;
    float tempR;
    float Massa[PDIM];
    float N, Media;
    float weight, simR;
    PLinkData linkD;


    N = GetNodes();
    RR = 0.0;
    
    //    if (verbose) cout << " Step#   \\theta_R   \\Delta\\theta_R \n";
    
    netDegree = 0.0;
    for (NI=BegNI(); NI<EndNI(); NI++) {
        data1 = GetNDat(NI.GetId());
        data1->degree = NI.GetOutDeg();
        netDegree += data1->degree;
        for (i=0 ; i<PDIM ; i++) {
            data1->xO[i] = data1->x[i];
            data1->dRO[i] = data1->dR[i];
        }
    }
    netDegree /= N;
    
    
    do {
        
        for (NI=BegNI(); NI<EndNI(); NI++) {
            data1 = GetNDat(NI.GetId());
            data1->DistF = 0;
            for (i=0 ; i<PDIM ; i++) {
                data1->dA[i] = 0;
                data1->dR[i] = 0;
                Massa[i] = 0.0;
            }
        }

        
        for (NI=BegNI() ; NI < EndNI() ; NI ++){
            id1 = NI.GetId();
            data1 = GetNDat(id1);
            d1 = data1->degree;
            for (NN=NI ; NN < EndNI() ; NN++){
                id2 = NN.GetId();
                if (id1==id2) continue;
                data2 = GetNDat(NN.GetId());
                d2 = data2->degree;

                
                for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - data2->x[i];
                r = 0;
                for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
                r = sqrt(r);
                simR = exp(-r);
                
                
                if (IsEdge(id1,id2)) {
                    linkD = GetEDat(id1,id2);
                    weight = linkD->weight;
                }
                else weight = 0.0;
                
                float cot;
                
                if (weight > simR) { // attraction
                    if (r<1.0) r=1.0; // to reduce the attraction when particles get close <1
                    for (i=0 ; i<PDIM ; i++) sum[i] = diff[i]/r;
                    for (i=0 ; i<PDIM ; i++) {
                        data1->dA[i] -= sum[i]*weight;
                        data2->dA[i] += sum[i]*weight;
                    }
                }
                else { // repulsion
                    data1->DistF += r;
                    data2->DistF += r;

                    if (r<0.001) r = 0.001; // to avoid division per zero when the particles are in the same position
                    for (i=0 ; i<PDIM ; i++) sum[i] = exp(-r)*diff[i]/(r);
                    
                    if (weight > 0) weight = 1.0 / weight;
                    else weight = 1000.0;
                    
                    weight = 1;
                    
                    for (i=0 ; i<PDIM ; i++) {
                        data1->dR[i] -= sum[i] * weight;
                        data2->dR[i] += sum[i] * weight;
                    }
                }
            }
        }
        
        
        oldRR = RR;
        RR = 0.0;
        for (NI=BegNI(); NI<EndNI(); NI++) {
            if ((degree = NI.GetOutDeg())){
                data1 = GetNDat(NI.GetId());
                tempR = 0;
                for (i=0 ; i<PDIM ; i++) {
                    //                    data1->xOL[i] = data1->x[i];
                    value = DT*((alpha*data1->dA[i] - beta*data1->dR[i]))/degree;
                    data1->x[i] += value;
                    data1->dAO[i] = value;
                    tempR += (data1->dR[i])*(data1->dR[i]);
                }
                RR += sqrt(tempR); //degree;
            }
        }
        // end Algorithm CORE
        
        ++steps;
        value = fabs(oldRR-RR); // / (float)Particles.size();
        if (verbose) cout << "Running: " << steps << " --> " << value << endl; //" DiffX " << diffX << endl;
//        if (value < minDR) {
//            if (steps%2) maxIT = steps;
//            else maxIT = steps+1;
//        }

        FILE * st;
        char str[256];
        
        sprintf(str,"Simulation.time.%d",steps);
        SaveParticlePosition(str);
    
    } while (steps<maxIT);
    
    
    
    
    // root-mean-square deviation of squared sums of repulsions and positions
//    diffX[0] = diffX[1] = diffX[2] = diffX[3] = 0.0;
//    for (NI=BegNI(); NI<EndNI(); NI++) {
//        data1 = GetNDat(NI.GetId());
//        degree = NI.GetOutDeg();
//        dTemp[0] = dTemp[1] = dTemp[2] = dTemp[3] = 0.0;
//        dTempR[0] = dTempR[1] = dTempR[2] = dTempR[3] = 0.0;
//        
//        for (i=0 ; i<PDIM ; i++) Massa[i] = 0.0;
//        
//        for (i=0 ; i<PDIM ; i++) {
//            Massa[i] += data1->x[i];
//            
//            value = data1->x[i] - data1->xO[i];
//            dTemp[0] += value*value;
//            dTemp[1] += value*value/degree;
//            dTemp[2] += value*value*degree;
//            if (degree > 1) dTemp[3] += value*value;
//            
//            value = data1->dR[i] - data1->dRO[i];
//            dTempR[0] += value*value;
//            dTempR[1] += value*value/degree;
//            dTempR[2] += value*value*degree;
//            if (degree>1)dTempR[3] += value*value;
//        }
//        diffX[0] += dTemp[0];
//        diffX[1] += dTemp[1];
//        diffX[2] += dTemp[2];
//        diffX[3] += dTemp[3];
//        
//        diffR[0] += dTempR[0];
//        diffR[1] += dTempR[1];
//        diffR[2] += dTempR[2];
//        diffR[3] += dTempR[3];
//        
//        for (i=0 ; i<PDIM ; i++) Massa[i] /= (N);
//        
//    }
//    DiffX[0] = diffX[0]/N;
//    DiffX[1] = diffX[1]/N;
//    DiffX[2] = diffX[2]/N;
//    DiffX[3] = diffX[3]/N;
//    
//    DiffR[0] = diffR[0]/N;
//    DiffR[1] = diffR[1]/N;
//    DiffR[2] = diffR[2]/N;
//    DiffR[3] = diffR[3]/N;
//    
//    
//    for (NI=BegNI(); NI<EndNI(); NI++) {
//        data1 = GetNDat(NI.GetId());
//        for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - Massa[i];
//        r = 0;
//        for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
//        r = sqrt(r);
//        data1->DistF = r; // distance from particle to the centre of mass
//        Media += r;
//    }
//    Media /= N;
//    
//    int ccd=0;
//    // mean-square deviation of distances from the centre of mass
//    for (NI=BegNI(); NI<EndNI(); NI++) {
//        data1 = GetNDat(NI.GetId());
//        degree = NI.GetOutDeg();
//        value = data1->DistF - Media;
//        VarMass[0] += value*value;
//        VarMass[1] += value*value/degree;
//        VarMass[2] += value*value*degree;
//        if (degree>1) {
//            VarMass[3] += value*value;
//            ccd++;
//        }
//    }
//    VarMass[0] /= N;
//    VarMass[1] /= N;
//    VarMass[2] /= N;
//    VarMass[3] /= (float)ccd;
    
    return steps;
}


int TParticleNet::sizeLargeCom(){
    int maxSize=0;
    vector<TCentroid>::iterator c, c_assigned;
    for (c=Centroids.begin() ; c!=Centroids.end(); ++c){
        if (c->nparticles > maxSize) maxSize = c->nparticles;
    }
    return maxSize;
}

void TParticleNet::CommunityDetectionDB(float epsilon){
    TParticleNet::TNodeI NI, NU;
    PParticle data1,data2;
    PLinkData linkD;
    int id1, id2;
    vector<int> neighbours;

    nextComDBId = 0;

    for (NI=BegNI(); NI<EndNI(); NI++) {
        data1 = GetNDat(NI.GetId());
        data1->visited = false;
        data1->cluster_db = 0;
        data1->neighbours = 0;
    }

    for (NI=BegNI(); NI<EndNI(); NI++) {
        id1 = NI.GetId();
        data1 = GetNDat(id1);
        if (data1->visited) continue;
        data1->visited = true;
        if (NI.GetDeg() == 1) neighbours.push_back(NI.GetOutNId(0));
        else
//        cout << "Debug1\n\n";
        for (i=0 ; i<NI.GetOutDeg() ; i++){
            id2 = NI.GetOutNId(i);
            data2 = GetNDat(id2);
            if (data2->visited) continue;
            linkD = GetEDat(id1,id2);
            if (linkD->distance < epsilon) neighbours.push_back(id2);
        }
        if (neighbours.size() > 2) {
//            cout << "Debug2\n\n";
            data1->cluster_db = ++nextComDBId;
            while (!neighbours.empty()){
                id1 = neighbours.back();
                neighbours.pop_back();
                data1 = GetNDat(id1);
                data1->visited = true;
                data1->cluster_db = nextComDBId;
                NU = GetNI(id1);

                for (i=0 ; i<NU.GetOutDeg() ; i++){
                    id2 = NU.GetOutNId(i);
                    data2 = GetNDat(id2);
                    if (data2->visited) continue;
                    linkD = GetEDat(id1,id2);
                    if (linkD->distance < epsilon) neighbours.push_back(id2);
                }
            }
        }
        else neighbours.clear();
    }
    


}




void TParticleNet::assignCentroids(){
    float dist, dist2, maxError=0.0, minError;
    vector<TCentroid>::iterator c, c_assigned;
    TParticleNet::TNodeI NI;
    PParticle p;
    int c_number=1;
    
    minError = numeric_limits<float>::max();

    for (c=Centroids.begin() ; c!=Centroids.end(); ++c){
        c->error = 0.0;
        c->nparticles = 0;
        c->comm_id = c_number++;
    }
    accError = 0.0;

    for (NI = BegNI() ; NI!=EndNI(); NI++){
        c = Centroids.begin();
        p = GetNDat(NI.GetId());

        dist = 0;
        for (i=0 ; i<PDIM ; i++) dist += pow(p->x[i] - c->x[i],2);

        c_assigned = c;
        ++c;
        while (c!=Centroids.end()){
            dist2 = 0;
            for (i=0 ; i<PDIM ; i++) dist2 += pow(p->x[i] - c->x[i],2);
            if (dist2 < dist){
                c_assigned = c;
                dist = dist2;
            }
            ++c;
        }

        p->index = &(*c_assigned);
        c_assigned->error += dist;
        c_assigned->nparticles++;
        if (dist > maxError) {
            idCentroidMaxError = c_assigned - Centroids.begin(); // Centroid's index associated to the farthest particle. When needed, the new centroid is added nearby this one.
            maxError = dist;
        }
        if (dist < minError){
            idCentroidMinError = c_assigned - Centroids.begin(); // Centroid's index associated to the farthest particle. When needed, the new centroid is added nearby this one.
            minError = dist;
            
        }
    }
    
    for (c=Centroids.begin() ; c!=Centroids.end(); ++c){
        accError += c->error;
    }
    

    toAddCentroid = false;
    toRemoveCentroid = false;
    for (c=Centroids.begin() ; c!=Centroids.end() ; ++c){
        if (c->nparticles) c->error /= c->nparticles;
        if (c->error == 0.0) {
            toRemoveCentroid = true;
            centroid2remove = c - Centroids.begin(); // centroid's index in the vector
        }
        else if (c->error > addCentroidThreshold){
            toAddCentroid = true;
        }
    }
    
}


void TParticleNet::computeCentroids(){
    TParticleNet::TNodeI NI;
    PParticle p;
    vector<TCentroid>::iterator c;
    
    for (c=Centroids.begin() ; c!=Centroids.end() ; ++c){
        for (i=0 ; i<PDIM ; i++) c->x[i] = 0;
        c->nparticles = 0;
        
        for (NI=BegNI() ; NI!=EndNI() ; NI++){
            p = GetNDat(NI.GetId());
            if (p->index == &(*c)){ //->comm_id){
                for (i=0 ; i<PDIM ; i++) c->x[i] += p->x[i];
                c->nparticles++;
            }
        }
        if (c->nparticles){
            for (i=0 ; i<PDIM ; i++) c->x[i] /= (float) c->nparticles;
        }
        else {
            for (i=0 ; i<PDIM ; i++) c->x[i] = 0.0;
        }
    }
}

               
void TParticleNet::addCentroid(){
    TCentroid centroid;
    for (i=0 ; i<PDIM ; i++) centroid.x[i] = Centroids[idCentroidMaxError].x[i] + (float)(rand()%10) / 100.0;
    centroid.comm_id = nextComId++;
    Centroids.push_back(centroid);
    toAddCentroid = false;
}

void TParticleNet::resetCentroid(){
    float vrand;
    for (i=0 ; i<PDIM ; i++) {
        vrand = -0.1 + (float)(rand()%201)/1000.0;
        Centroids[idCentroidMinError].x[i] = Centroids[idCentroidMaxError].x[i] + vrand;
    }
    
}

void TParticleNet::removeCentroid(){
    Centroids.erase(Centroids.begin()+centroid2remove);
    toRemoveCentroid = false;
}

void TParticleNet::mergeCentroids(){
    vector<TCentroid>::iterator c1,c2,c_remove;
    float min=99999, dist;

    centroidsMerged = false;
    for (c1=Centroids.begin() ; c1!=(Centroids.end()-1) ; ++c1){
        for (c2=c1+1 ; c2!=(Centroids.end()) ; ++c2){
            dist = 0;
            for (i=0 ; i<PDIM ; i++) dist += pow(c1->x[i] - c2->x[i],2);
            if (dist < min){
                min = dist;
                c_remove = c1;
            }
        }
    }
    if (sqrt(min)<0.9) {
        centroidsMerged = true;
        Centroids.erase(c_remove);
//        cout << "M";
    }
}

int TParticleNet::CommunityDetection3(){
    float oldAcc;
    float diffError;
    int steps=0, countDown=5000;
    
    toAddCentroid = false;
    toRemoveCentroid = false;
    
    do {
//        cout << "Step: " << steps << " ";
        oldAcc = accError;
        mergeCentroids();
        assignCentroids();
        if (toAddCentroid || toRemoveCentroid){
            if (toRemoveCentroid) {
                removeCentroid();
//                cout << "R";
            }
            if (toAddCentroid) {
                addCentroid();
//                cout << "A";
            }
            assignCentroids();
        }
        computeCentroids();
        accError = accError*0.1 + oldAcc*0.9; // used to make the evolution of the centroid error more stable...
        diffError = fabs(oldAcc - accError); // / (float)Centroids.size();
        ++steps;
        countDown--;
//        cout << " #:" << Centroids.size() << endl;
    } while (diffError>0.01 && countDown);
    return steps;
}


float TParticleNet::printCentroidsError(){
    return accError;
}


//bool myfunction (TParticle i, TParticle j) { return (i.node_id<j.node_id); }
//
//void TParticleNet::SaveCentroids(const char *filename){
//    ofstream file;
//        
//    vector<TCentroid>::iterator i;
//
//    int aux;
//    for (aux=1, i=Centroids.begin() ; i!=Centroids.end(); ++i, aux++){
//        i->comm_id = aux;
//    }
//    
//    file.open(filename, ofstream::out);
//    file << "% id_centroids; position; # of associated particles; total repulsion; total attraction; sum of dists " << endl;
//    for (i = Centroids.begin() ; i != Centroids.end(); ++i){
//        file << i->x << "\t"
//             << i->y << "\t"
//             << i->z << "\t"
//             << (int) i->comm_id << "\t"
//             << i->nparticles << "\t"
//             << i->totalR << "\t"
//             << i->totalA << "\t"
//             << i->error << endl;
//    }
//    file.close();
//}
//
//    
//void TParticleNet::SaveCommunities(const char *filename){
//    ofstream file;
//    
//    int aux;
//    vector<TCentroid>::iterator c;
//    for (aux=1, c=Centroids.begin() ; c!=Centroids.end(); ++c, aux++){
//        c->comm_id = aux;
//    }
//    sort (Particles.begin(), Particles.end(), myfunction);
//    
//    vector<TParticle>::iterator i;
//    
//    file.open(filename, ofstream::out);
//    for (i = Particles.begin() ; i != Particles.end(); ++i){
////        file << i->index << "\t" << i->node_id << endl;
//        file << i->node_id <<  "\t" << i->index->comm_id << endl;
//    }
//    file.close();
//}
//
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


float TParticleNet::NMI(){
    int **labelCount;
    int *numCount, *sumI, *sumJ;
    float precision, Si=0,Sj=0;
    int p, valor, comUser;
    int N, numCentroids;
    
    if (Centroids.size()<=1) return 0;
    
    vector<TCentroid>::iterator c;

    for (i=0, c=Centroids.begin() ; c!=Centroids.end(); ++c, i++){
        c->comm_id = i;
    }
    
//    	cout << "Com: " << numCommunities << endl;
//    	cout << "Cen: " << numCentroids << endl;
    
    numCentroids = Centroids.size();
    N = GetNodes();
    
    if (numCentroids > numCommunities) comUser = numCentroids;
    else comUser = numCommunities;
    
    labelCount = new int*[comUser];
    sumI = new int[comUser];
    sumJ = new int[comUser];
    for (i=0 ; i<comUser ; i++){
        sumI[i] = sumJ[i] = 0.0;
        labelCount[i] = new int[comUser];
        for (j=0 ; j<comUser ; j++){
            labelCount[i][j] = 0;
        }
    }

    TParticleNet::TNodeI NI;
    PParticle particle;
    
    for (NI=BegNI(); NI<EndNI(); NI++) {
        particle = GetNDat(NI.GetId());
        particle->refreshed = false;
        labelCount[particle->index->comm_id][particle->indexReal-1]++;
//        cout << NI.GetId() << " --> " << NI.GetInDeg() << endl;

    }

    
    // sorting by rows
    for (j=0 ; j<comUser-1 ; j++){
        p = j;
        valor = labelCount[j][j];
        for (i=j+1 ; i<comUser ; i++){
            if (labelCount[i][j] > valor){
                p = i;
                valor = labelCount[i][j];
            }
        }
        if (p!=j) {
            numCount = labelCount[j];
            labelCount[j] = labelCount[p];
            labelCount[p] = numCount;
        }
    }
    
    // sorting by columns
    for (j=0 ; j<comUser-1 ; j++){
        p = j;
        valor = labelCount[j][j];
        for (i=j+1 ; i<comUser ; i++){
            if (labelCount[j][i] > valor){
                p = i;
                valor = labelCount[j][i];
            }
        }
        if (p!=j) {
            for (i=0 ; i<comUser ; i++){
                valor = labelCount[i][j];
                labelCount[i][j] = labelCount[i][p];
                labelCount[i][p] = valor;
            }
        }
    }
    

//    cout << "Matrix:\n";
    for (j=0 ; j<comUser ; j++){
        for (i=0 ; i<comUser ; i++){
            sumJ[j] += labelCount[j][i];
            sumI[i] += labelCount[j][i];
//            cout << labelCount[j][i] << "\t";
        }
//        cout << ";" << endl;
    }

    comUser = numCommunities;
    
    precision = 0.0;
    Si = 0.0;
    Sj = 0.0;
    int pNormal=0;
    for (i=0 ; i<comUser ; i++){
        pNormal += labelCount[i][i];
        for (j=0 ; j<comUser ; j++){
            if (labelCount[j][i]){
                precision += (float)labelCount[j][i] * log( (float)(labelCount[j][i]*N)  /  (float)(sumI[i]*sumJ[j]) );
            }
        }
        if (sumI[i]) Si += (float)sumI[i]*log((float)sumI[i]/(float)N);
        if (sumJ[i]) Sj += (float)sumJ[i]*log((float)sumJ[i]/(float)N);
    }
    
    precision *= -2;
    precision /= (Si+Sj);
    
//    for (i=0 ; i<comUser ; i++)
//        delete[] labelCount[i];
//    delete[] labelCount;
//    delete[] sumI;
//    delete[] sumJ;
    
    return precision;
}

//float TParticleNet::NMI2(){
//    vector<TParticle>::iterator part;
//    int i,j;
//    int **labelCount;
//    int *numCount, *sumI, *sumJ;
//    float precision, Si=0,Sj=0;
//    int p, valor, comUser;
//    int N;
//    
//    if (!numClusters) return 0;
//    
////    numClusters--;
//    
//    N = Particles.size();
//
//    if (numClusters > numCommunities) comUser = numClusters;
//    else comUser = numCommunities;
//    
////    cout << "#Clusters: " << numClusters << endl;
//
//    labelCount = new int*[comUser];
//    sumI = new int[comUser];
//    sumJ = new int[comUser];
//    for (i=0 ; i<comUser ; i++){
//        sumI[i] = sumJ[i] = 0.0;
//        labelCount[i] = new int[comUser];
//        for (j=0 ; j<comUser ; j++){
//            labelCount[i][j] = 0;
//        }
//    }
//
////    cout << "#Clusters: " << numClusters << endl;
//    
//    
//    for (part = Particles.begin() ; part != Particles.end() ; ++part){
////        cout << part->node_id << " -> " << part->cluster_id << endl;
////        labelCount[part->cluster_id-1][part->indexReal-1]++;
//        labelCount[part->cluster_id][part->indexReal-1]++;
//    }
//    
//    
//    // sorting by rows
//    for (j=0 ; j<comUser-1 ; j++){
//        p = j;
//        valor = labelCount[j][j];
//        for (i=j+1 ; i<comUser ; i++){
//            if (labelCount[i][j] > valor){
//                p = i;
//                valor = labelCount[i][j];
//            }
//        }
//        if (p!=j) {
//            numCount = labelCount[j];
//            labelCount[j] = labelCount[p];
//            labelCount[p] = numCount;
//        }
//    }
//    
//    // sorting by columns
//    for (j=0 ; j<comUser-1 ; j++){
//        p = j;
//        valor = labelCount[j][j];
//        for (i=j+1 ; i<comUser ; i++){
//            if (labelCount[j][i] > valor){
//                p = i;
//                valor = labelCount[j][i];
//            }
//        }
//        if (p!=j) {
//            for (i=0 ; i<comUser ; i++){
//                valor = labelCount[i][j];
//                labelCount[i][j] = labelCount[i][p];
//                labelCount[i][p] = valor;
//            }
//        }
//    }
//    
//    
////    cout << "Matrix:\n";
//    for (j=0 ; j<comUser ; j++){
//        for (i=0 ; i<comUser ; i++){
//            sumJ[j] += labelCount[j][i];
//            sumI[i] += labelCount[j][i];
////            cout << labelCount[j][i] << "\t";
//        }
////        cout << ";" << endl;
//    }
//    
//    comUser = numCommunities;
//    
//    precision = 0.0;
//    Si = 0.0;
//    Sj = 0.0;
//    int pNormal=0;
//    for (i=0 ; i<comUser ; i++){
//        pNormal += labelCount[i][i];
//        for (j=0 ; j<comUser ; j++){
//            if (labelCount[j][i]){
//                precision += (float)labelCount[j][i] * log( (float)(labelCount[j][i]*N)  /  (float)(sumI[i]*sumJ[j]) );
//            }
//        }
//        if (sumI[i]) Si += (float)sumI[i]*log((float)sumI[i]/(float)N);
//        if (sumJ[i]) Sj += (float)sumJ[i]*log((float)sumJ[i]/(float)N);
//    }
//    
//    precision *= -2;
//    precision /= (Si+Sj);
//    
//    //    for (i=0 ; i<comUser ; i++)
//    //        delete[] labelCount[i];
//    //    delete[] labelCount;
//    //    delete[] sumI;
//    //    delete[] sumJ;
//    
//    return precision;
//}
//
//float TParticleNet::NMIH(int nivel){
//    
//    int i,j;
//    int **labelCount;
//    int *numCount, *sumI, *sumJ;
//    float precision, Si=0,Sj=0;
//    int p, valor, comUser;
//    int N, numCentroids;
//    
//    
//    ShrinkCentroids();
//    
//    if (Centroids.size()<=1) return 0;
//    
//    vector<TCentroid>::iterator c;
//    vector<TParticle>::iterator part;
//    for (i=0, c=Centroids.begin() ; c!=Centroids.end(); ++c, i++){
//        c->comm_id = i;
//    }
//    
//    //    	cout << "Com: " << numCommunities << endl;
//    //    	cout << "Cen: " << numCentroids << endl;
//    
//    numCentroids = Centroids.size();
//    N = Particles.size();
//    
//    numCommunities = numCommunitiesH[nivel];
//    
////    cout << "\n\n# Com: " << numCommunities << endl;
//    
//    if (numCentroids > numCommunities) comUser = numCentroids;
//    else comUser = numCommunities;
//    
//    labelCount = new int*[comUser];
//    sumI = new int[comUser];
//    sumJ = new int[comUser];
//    for (i=0 ; i<comUser ; i++){
//        sumI[i] = sumJ[i] = 0.0;
//        labelCount[i] = new int[comUser];
//        for (j=0 ; j<comUser ; j++){
//            labelCount[i][j] = 0;
//        }
//    }
//    
//    for (part = Particles.begin() ; part != Particles.end() ; ++part){
//        labelCount[part->index->comm_id][part->indexRealH[nivel]-1]++;
////        cout << "Par: " << part->node_id << " com ac: " << part->index->comm_id << " com re: " << part->indexRealH[nivel] << endl;
//    }
//    
//    
//    // sorting by rows
//    for (j=0 ; j<comUser-1 ; j++){
//        p = j;
//        valor = labelCount[j][j];
//        for (i=j+1 ; i<comUser ; i++){
//            if (labelCount[i][j] > valor){
//                p = i;
//                valor = labelCount[i][j];
//            }
//        }
//        if (p!=j) {
//            numCount = labelCount[j];
//            labelCount[j] = labelCount[p];
//            labelCount[p] = numCount;
//        }
//    }
//    
//    // sorting by columns
//    for (j=0 ; j<comUser-1 ; j++){
//        p = j;
//        valor = labelCount[j][j];
//        for (i=j+1 ; i<comUser ; i++){
//            if (labelCount[j][i] > valor){
//                p = i;
//                valor = labelCount[j][i];
//            }
//        }
//        if (p!=j) {
//            for (i=0 ; i<comUser ; i++){
//                valor = labelCount[i][j];
//                labelCount[i][j] = labelCount[i][p];
//                labelCount[i][p] = valor;
//            }
//        }
//    }
//    
//    
////        cout << "Matrix:\n";
//    for (j=0 ; j<comUser ; j++){
//        for (i=0 ; i<comUser ; i++){
//            sumJ[j] += labelCount[j][i];
//            sumI[i] += labelCount[j][i];
////                        cout << labelCount[j][i] << "\t";
//        }
////                cout << ";" << endl;
//    }
////    cout << endl;
//    
//    comUser = numCommunities;
//    
//    precision = 0.0;
//    Si = 0.0;
//    Sj = 0.0;
//    int pNormal=0;
//    for (i=0 ; i<comUser ; i++){
//        pNormal += labelCount[i][i];
//        for (j=0 ; j<comUser ; j++){
//            if (labelCount[j][i]){
//                precision += (float)labelCount[j][i] * log( (float)(labelCount[j][i]*N)  /  (float)(sumI[i]*sumJ[j]) );
//            }
//        }
//        if (sumI[i]) Si += (float)sumI[i]*log((float)sumI[i]/(float)N);
//        if (sumJ[i]) Sj += (float)sumJ[i]*log((float)sumJ[i]/(float)N);
//    }
//    
//    precision *= -2;
//    precision /= (Si+Sj);
//    
//    //    for (i=0 ; i<comUser ; i++)
//    //        delete[] labelCount[i];
//    //    delete[] labelCount;
//    //    delete[] sumI;
//    //    delete[] sumJ;
//    
//    return precision;
//}
//
//
//
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
//
//
//// fixed << setprecision(2)
//

float TParticleNet::getNormFR(int d){
    TParticleNet::TNodeI NI;
    PParticle data;
    float DR, Total=0.0, value;
    float degree;
    
    for (NI=BegNI() ; NI<EndNI(); NI++){
        data = GetNDat(NI.GetId());
        degree = NI.GetOutDeg();
        DR = 0.0;
        for (i=0 ; i<PDIM ; i++) {
            value = data->dR[i];
            if (d==0) DR += value*value;
            else if (d==1) DR += value*value/degree;
            else if (d==2) DR += value*value*degree;
            else if (d==3) {
                if (degree>1) DR += value*value;
            }
        }
//        DR = sqrt(DR);
        Total += DR; //pow(DR,2);//DR*DR;
    }
    return Total;
}


/*
 PUNGraph UGraph = TSnap::ConvertGraph<PUNGraph>(Graph); // undirected version of the graph
 TIntFltH BtwH, EigH, PRankH, CcfH, CloseH, HubH, AuthH;
 //printf("Computing...\n");
 printf("Treat graph as DIRECTED: ");
 printf(" PageRank... ");             TSnap::GetPageRank(Graph, PRankH, 0.85);
 printf(" Hubs&Authorities...");      TSnap::GetHits(Graph, HubH, AuthH);
 printf("\nTreat graph as UNDIRECTED: ");
 printf(" Eigenvector...");           TSnap::GetEigenVectorCentr(UGraph, EigH);
 printf(" Clustering...");            TSnap::GetNodeClustCf(UGraph, CcfH);
 printf(" Betweenness (SLOW!)...");   TSnap::GetBetweennessCentr(UGraph, BtwH, 1.0);
 printf(" Constraint (SLOW!)...");    TNetConstraint<PUNGraph> NetC(UGraph, true);
 printf(" Closeness (SLOW!)...");
 for (TUNGraph::TNodeI NI = UGraph->BegNI(); NI < UGraph->EndNI(); NI++) {
 const int NId = NI.GetId();
 CloseH.AddDat(NId, TSnap::GetClosenessCentr<PUNGraph>(UGraph, NId, false));
 }
 printf("\nDONE! saving...");
 FILE *F = fopen(OutFNm.CStr(), "wt");
 fprintf(F,"#Network: %s\n", InFNm.CStr());
 fprintf(F,"#Nodes: %d\tEdges: %d\n", Graph->GetNodes(), Graph->GetEdges());
 fprintf(F,"#NodeId\tDegree\tCloseness\tBetweennes\tEigenVector\tNetworkConstraint\tClusteringCoefficient\tPageRank\tHubScore\tAuthorityScore\n");
 for (TUNGraph::TNodeI NI = UGraph->BegNI(); NI < UGraph->EndNI(); NI++) {
 const int NId = NI.GetId();
 const double DegCentr = UGraph->GetNI(NId).GetDeg();
 const double CloCentr = CloseH.GetDat(NId);
 const double BtwCentr = BtwH.GetDat(NId);
 const double EigCentr = EigH.GetDat(NId);
 const double Constraint = NetC.GetNodeC(NId);
 const double ClustCf = CcfH.GetDat(NId);
 const double PgrCentr = PRankH.GetDat(NId);
 const double HubCentr = HubH.GetDat(NId);
 const double AuthCentr = AuthH.GetDat(NId);
 fprintf(F, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", NId,
 DegCentr, CloCentr, BtwCentr, EigCentr, Constraint, ClustCf, PgrCentr, HubCentr, AuthCentr);
 }

 */

void TParticleNet::SaveMeasures(const char *filename){
//    ofstream file;
//    TParticleNet::TNodeI NI;
//    
//    PParticleNet graph;
//    graph = this;
////
////    //PUNGraph UGraph = TSnap::ConvertGraph<PUNGraph>(Graph); // undirected version of the graph
//    TIntFltH BtwH, EigH, PRankH, CcfH, CloseH, HubH, AuthH;
//    //printf("Computing...\n");
////    printf("Treat graph as DIRECTED: ");
////    printf(" PageRank... ");             TSnap::GetPageRank(Graph, PRankH, 0.85);
////    printf(" Hubs&Authorities...");      TSnap::GetHits(Graph, HubH, AuthH);
////    printf("\nTreat graph as UNDIRECTED: ");
//    printf(" Eigenvector...");           TSnap::GetEigenVectorCentr(graph, EigH);
//    printf(" Clustering...");            TSnap::GetNodeClustCf(this, CcfH);
////    printf(" Betweenness (SLOW!)...");   TSnap::GetBetweennessCentr(UGraph, BtwH, 1.0);
////    printf(" Constraint (SLOW!)...");    TNetConstraint<PUNGraph> NetC(UGraph, true);
////    printf(" Closeness (SLOW!)...");
//
//    
//    file.open(filename, ofstream::out);
//    file << "% Particle's file -> a snapshot of the particle space\n";
//    file << "% Simulation Parameters -> alpha: " << alpha << " beta: " << beta << "Particle dimension: " << PDIM << endl;
//    file << "% node_id, ground truth community id, assigned community id, node_id, x1, x2, fa(x1), fa(x2), ..., fr(x1), fr(x2), ...\n";
//    
//    for (NI=BegNI() ; NI<EndNI(); NI++){
//        
////            file << fixed << setprecision(3) << DA << "\t" << DR << "\t";
////            file << endl;
//    }
//    file.close();
}

void TParticleNet::SaveParticleForces(const char *filename){
    ofstream file;
    TParticleNet::TNodeI NI;
    PParticle data;
    float DR, DA, DAO, DRO; 

    file.open(filename, ofstream::out);
    file << "% Particle's file -> a snapshot of the particle space\n";
    file << "% Simulation Parameters -> alpha: " << alpha << " beta: " << beta << " Particle dimension: " << PDIM << endl;
    file << "% node_id, degree, dist(t=0), dist(t=inf), Atraction (t=0), Atraction (t=inf), Repulsion (t=0), Repulsion (t=inf)\n";

    for (NI=BegNI() ; NI<EndNI(); NI++){
        data = GetNDat(NI.GetId());
        
        file << fixed << setprecision(3) << NI.GetId() << "\t" << NI.GetOutDeg() << "\t" << data->DistI << "\t" << data->DistF << "\t";
        
        DA = DR = DAO = DRO = 0.0;
        for (i=0 ; i<PDIM ; i++) {
            DA += pow(data->dA[i],2.0);
            DR += pow(data->dR[i],2.0);
            DAO += pow(data->dAO[i],2.0);
            DRO += pow(data->dRO[i],2.0);
        }
        
        file << fixed << setprecision(3) << sqrt(DAO) << "\t" << sqrt(DA)  << "\t" << sqrt(DRO)  << "\t" << sqrt(DR) << endl;
        
    }
    file.close();
}


void TParticleNet::SaveParticlePosition(const char *filename){
    ofstream file;
    TParticleNet::TNodeI NI;
    PParticle data;
    float DR, DA, DAO, DRO; 

    file.open(filename, ofstream::out);
    file << "% Particle's file -> a snapshot of the particle space\n";
    file << "% Simulation Parameters -> alpha: " << alpha << " beta: " << beta << " Particle dimension: " << PDIM << endl;
//    file << "% node_id, node degree, centrality, ground truth community id, assigned com id1, assigned com dbscan x1, x2, ..., x_dim\n";

//    file << "@RELATION Iris\n\n";
//    file << "@ATTRIBUTE d1 REAL\n";
//    file << "@ATTRIBUTE d2 REAL\n";
//    file << "@ATTRIBUTE d3 REAL\n\n@DATA\n";
    
    for (NI=BegNI() ; NI<EndNI(); NI++){
        data = GetNDat(NI.GetId());
        
        if (data->index) {
            file << fixed << setprecision(2) << NI.GetId() << "\t " // 1
                                             << NI.GetOutDeg() << "\t" //2
//                                             << data->ComCentrality << "\t" //3
                                             << data->indexReal << "\t"  //4
                                             << data->index->comm_id << "\t"; //5
//                                             << data->cluster_db << "\t"; //6
        }
        else {
            file << fixed << setprecision(2) << NI.GetId() << "\t " 
                                             << NI.GetOutDeg() << "\t"
//                                             << data->ComCentrality << "\t"
                                             << data->indexReal << "\t"
                                             << "-1" << "\t";
//                                             << data->cluster_db << "\t";
        }
//        file << fixed << setprecision(2) << data->x[0] << ", ";
//        file << fixed << setprecision(2) << data->x[1] << ", ";
//        file << fixed << setprecision(2) << data->x[2] << "\n";
        for (i=0 ; i<PDIM ; i++)  file << fixed << setprecision(2) << data->x[i] << "\t";
//        for (i=0 ; i<PDIM ; i++)  file << fixed << setprecision(2) << data->xO[i] << "\t";
//        for (i=0 ; i<PDIM ; i++)  file << fixed << setprecision(2) << data->dA[i] << "\t";
//        for (i=0 ; i<PDIM ; i++)  file << fixed << setprecision(2) << data->dR[i] << "\t";
//        for (i=0 ; i<PDIM ; i++)  file << fixed << setprecision(2) << data->dAO[i] << "\t";
        file << "\n";
    }
    file.close();
}
//
//int TParticleNet::getNumCommunities(){
//    return Centroids.size();
//}
//
//int TParticleNet::getNumParticles(){
//    return Particles.size();
//}
//
//void TParticleNet::fineTuning(int t){
//    int i;
//    for (i=0 ; i<t ; i++){
//        assignCentroids();
//        computeCentroids();
//    }
//}
//
//void TParticleNet::ShrinkCentroids(){
//    vector<TCentroid>::iterator i;
//    int aux;
//    for (aux=1, i=Centroids.begin() ; i!=Centroids.end(); ++i, aux++){
//        i->comm_id = aux;
//    }
//}
//
//void TParticleNet::printCentroids(){
//    vector<TCentroid>::iterator i;
//    float E=0.0, TR=0.0, TA=0.0;
//    
//    int aux;
//    for (aux=1, i=Centroids.begin() ; i!=Centroids.end(); ++i, aux++){
//        i->comm_id = aux;
//    }
//    
//    for (i = Centroids.begin() ; i != Centroids.end(); ++i){
////        cout << "C: " << i->comm_id << " E: " << i->error << " A: " << i->totalA << " R: " << i->totalR << endl;
//        E += i->error;
//        TA += i->totalA;
//        TR += i->totalR;
////        file << i->x << "\t"
////        << i->y << "\t"
////        << i->z << "\t"
////        << (int) i->comm_id << "\t"
////        << i->nparticles << "\t"
////        << i->totalR << "\t"
////        << i->totalA << "\t"
////        << i->error << endl;
//    }
////    file.close();
//    E /= Centroids.size();
//    TA /= Centroids.size();
//    TR /= Centroids.size();
//    cout << " " << E << " " << TA << " " << TR << endl;
////    cout << " " << E << " " << Centroids.size() << endl;
//
//}
//
//float TParticleNet::printCentroidsError(){
//    return accError;
//}
//
//
//////////////////
//////////////////
//////////////////
//////////////////
//////////////////
//////////////////
//////////////////
//////////////////
//////////////////
//////////////////
//////////////////
//
//
//#define EPS 1.0
//#define MINPTS 4
//
//vector<TParticle*> TParticleNet::GetRegion(TParticle *p){
//    float dist;
//    vector<TParticle>::iterator n;
//    vector<TParticle*> region;
//    
//    for (n=Particles.begin() ; n!=Particles.end(); ++n){
//        if (&(*n)==p) continue;
//        dist = pow(n->x - p->x,2) + pow(n->y - p->y,2) + pow(n->z - p->z,2);
//        dist = sqrt(dist);
//        if (dist <= EPS) region.push_back(&(*n));
////        cout << dist << "\t";
//    }
//    return region;
//}
//
//bool TParticleNet::ExpandCluster(TParticle *p, int clusterID){
//    
//    vector<TParticle*> seeds = GetRegion(p);
//    
//    if (seeds.size() < MINPTS)
//    {
//        p->cluster_id = -1;
//        return false;
//    }
//    else // all points in seeds are density reachable from point 'p'
//    {
//        p->cluster_id = clusterID;
//        for (int s=0; s<seeds.size(); s++) {
//            seeds[s]->cluster_id = clusterID;
////            cout << " " << seeds[s]->node_id;
//        }
//        while (seeds.size() > 0)
//        {
//            TParticle* currentP = seeds.back();
//            seeds.pop_back();
//            vector<TParticle*> seeds2 = GetRegion(currentP);
//            if (seeds2.size() >= MINPTS)
//            {
//                for (int i=0; i<seeds2.size(); i++)
//                {
//                    TParticle* currentP2 = seeds2.back();
//                    if (currentP2->cluster_id == 0  || currentP2->cluster_id == -1)
//                    {
//                        if (currentP2->cluster_id == 0) seeds.push_back(currentP2);
//                        currentP2->cluster_id = clusterID;
////                        cout << " " << currentP2->node_id;
//                    }
//                }
//            }
//        }
//        return true;
//    }
//}
//
//void TParticleNet::CommunityDetection2(){
//    
//    vector<TParticle>::iterator i;
//    int clusterID=1;
//
//    for (i=Particles.begin() ; i!=Particles.end(); ++i){
//        i->cluster_id = 0;
//    }
//    for (i = Particles.begin() ; i != Particles.end(); ++i){
////        cout << "Seed: " << i->node_id;
//        if (i->cluster_id == 0){
//            if (ExpandCluster(&(*i),clusterID)) clusterID++;
//        }
////        cout << endl;
//    }
//    numClusters = clusterID;
//
////    for (i=Particles.begin() ; i!=Particles.end(); ++i){
////        cout << "Particle: " << i->cluster_id << endl;
////    }
//    for (i=Particles.begreset clockin() ; i!=Particles.end(); ++i){
//        if (i->cluster_id == -1) i->cluster_id = 0;
//    }
//
//
//    
//}



//void TParticleNet::AddNoise(){
//    float sumX, sumY, sumZ, r;
//    float diffx, diffy, diffz;
//    TParticleNet::TNodeI NI, NN;
//    TParticleNet::TEdgeI EI;
//    PParticle data1, data2;
//    int id1, id2, degree;
//
//    sumX = sumY = sumZ = 0.0;
//    for (NI=BegNI(); NI<EndNI(); NI++) {
//        if ((degree = NI.GetDeg())){
//            data1 = GetNDat(NI.GetId());
//            data1->x += DT*((alpha*data1->dxA - beta*data1->dxR)) / (float)degree;
//            data1->y += DT*((alpha*data1->dyA - beta*data1->dyR)) / (float)degree;
//            data1->z += DT*((alpha*data1->dzA - beta*data1->dzR)) / (float)degree;
//            
//            // Debug
//            data1->dxA = (alpha*data1->dxA)/(float)degree;
//            data1->dyA = (alpha*data1->dyA)/(float)degree;
//            data1->dzA = (alpha*data1->dzA)/(float)degree;
//            data1->dxR = (beta*data1->dxR)/(float)degree;
//            data1->dyR = (beta*data1->dyR)/(float)degree;
//            data1->dzR = (beta*data1->dzR)/(float)degree;
//        }
//    }
//    
//}

/*
void TParticleNet::RunByStepRadial(){
    float sumX, sumY, sumZ, r;
    float diffx, diffy, diffz;
    TParticleNet::TNodeI NI, NN;
    TParticleNet::TEdgeI EI;
    PParticle data1, data2;
    int id1, id2, degree;
    
    for (NI=BegNI(); NI<EndNI(); NI++) {
        data1 = GetNDat(NI.GetId());
        data1->dxA = 0;
        data1->dyA = 0;
        data1->dzA = 0;
        data1->dxR = 0;
        data1->dyR = 0;
        data1->dzR = 0;
    }


    // Attraction
    for (EI=BegEI() ; EI < EndEI() ; EI++){
        id1 = EI.GetSrcNId();
        id2 = EI.GetDstNId();
        // Loading node's data;
        data1 = GetNDat(id1);
        data2 = GetNDat(id2);
        diffx = data1->x - data2->x;
        diffy = data1->y - data2->y;
        diffz = data1->z - data2->z;
        r = pow(diffx,2) + pow(diffy,2) + pow(diffz,2);
        r = sqrt(r);
        if (r<1.0) r = 4.0;
        sumX = diffx/r;
        sumY = diffy/r;
        sumZ = diffz/r;
        data1->dxA -= sumX;
        data1->dyA -= sumY;
        data1->dzA -= sumZ;
        data2->dxA += sumX;
        data2->dyA += sumY;
        data2->dzA += sumZ;
    }

    float massX=0, massY=0, massZ=0;
    for (NI=BegNI() ; NI < EndNI() ; NI ++){
        data1 = GetNDat(NI.GetId());
        massX += data1->x;
        massY += data1->y;
        massZ += data1->z;
    }
    massX /= 128.0;
    massY /= 128.0;
    massZ /= 128.0;
    
    //Repulsion
    for (NI=BegNI() ; NI < EndNI() ; NI ++){
        data1 = GetNDat(NI.GetId());

        diffx = data1->x - massX;
        diffy = data1->y - massY;
        diffz = data1->z - massZ;
        r = pow(diffx,2) + pow(diffy,2) + pow(diffz,2);
        r = sqrt(r);

        sumX = exp(-1.0*(r-1.0))*diffx/(r);
        sumY = exp(-1.0*(r-1.0))*diffy/(r);
        sumZ = exp(-1.0*(r-1.0))*diffz/(r);
        data1->dxR += sumX;
        data1->dyR += sumY;
        data1->dzR += sumZ;
//            data2->dxR += sumX;
//            data2->dyR += sumY;
//            data2->dzR += sumZ;

    }
    
//    alpha = 1.0;
//    beta = 0.2;
    sumX = sumY = sumZ = 0.0;
    for (NI=BegNI(); NI<EndNI(); NI++) {
        if ((degree = NI.GetDeg())){
            data1 = GetNDat(NI.GetId());
            data1->x += ((alpha*data1->dxA - beta*data1->dxR)) / (float)degree;
            data1->y += ((alpha*data1->dyA - beta*data1->dyR)) / (float)degree;
            data1->z += ((alpha*data1->dzA - beta*data1->dzR)) / (float)degree;
        }
    }
}

void TParticleNet::RunByStepRadial2(){
    float sumX, sumY, sumZ, r;
    float diffx, diffy, diffz;
    TParticleNet::TNodeI NI, NN;
    TParticleNet::TEdgeI EI;
    PParticle data1, data2;
    int id1, id2, degree;
    
    for (NI=BegNI(); NI<EndNI(); NI++) {
        data1 = GetNDat(NI.GetId());
        data1->dxA = 0;
        data1->dyA = 0;
        data1->dzA = 0;
        data1->dxR = 0;
        data1->dyR = 0;
        data1->dzR = 0;
        data1->AxA = 0;
        data1->AyA = 0;
        data1->AzA = 0;
        data1->AxR = 0;
        data1->AyR = 0;
        data1->AzR = 0;

    }
    
    
    // Attraction
    for (EI=BegEI() ; EI < EndEI() ; EI++){
        id1 = EI.GetSrcNId();
        id2 = EI.GetDstNId();
        // Loading node's data;
        data1 = GetNDat(id1);
        data2 = GetNDat(id2);
        diffx = data1->x - data2->x;
        diffy = data1->y - data2->y;
        diffz = data1->z - data2->z;
        r = pow(diffx,2) + pow(diffy,2) + pow(diffz,2);
        r = sqrt(r);
//        r = 1.0 / (1.0 + exp(-r));
//        if (r==0) r = 0.00001; // to avoid division per zero when two particles are in the same position
//        r = 1.0;
//        sumX = diffx*r;
//        sumY = diffy*r;
//        sumZ = diffz*r;
//        sumX = diffx/r;
//        sumY = diffy/r;
//        sumZ = diffz/r;
        sumX = (diffx/r)*(1.0 + 0.0/r);
        sumY = (diffy/r)*(1.0 + 0.0/r);
        sumZ = (diffz/r)*(1.0 + 0.0/r);
        data1->dxA -= sumX;
        data1->dyA -= sumY;
        data1->dzA -= sumZ;
        data2->dxA += sumX;
        data2->dyA += sumY;
        data2->dzA += sumZ;

        data1->AxA -= fabs(sumX);
        data1->AyA -= fabs(sumY);
        data1->AzA -= fabs(sumZ);
        data2->AxA += fabs(sumX);
        data2->AyA += fabs(sumY);
        data2->AzA += fabs(sumZ);
    }
    
    float massX=0, massY=0, massZ=0;
    for (NI=BegNI() ; NI < EndNI() ; NI ++){
        data1 = GetNDat(NI.GetId());
        massX += data1->x;
        massY += data1->y;
        massZ += data1->z;
    }
    massX /= 128.0;
    massY /= 128.0;
    massZ /= 128.0;
    
    //Repulsion
    for (NI=BegNI() ; NI < EndNI() ; NI ++){
        data1 = GetNDat(NI.GetId());
        
        diffx = data1->x - massX;
        diffy = data1->y - massY;
        diffz = data1->z - massZ;
        r = pow(diffx,2) + pow(diffy,2) + pow(diffz,2);
//        r = sqrt(r);
        sumX = diffx/(r);
        sumY = diffy/(r);
        sumZ = diffz/(r);
        data1->dxR -= sumX;
        data1->dyR -= sumY;
        data1->dzR -= sumZ;
        data2->dxR += sumX;
        data2->dyR += sumY;
        data2->dzR += sumZ;
        
        data1->AxR -= fabs(sumX);
        data1->AyR -= fabs(sumY);
        data1->AzR -= fabs(sumZ);
        data2->AxR += fabs(sumX);
        data2->AyR += fabs(sumY);
        data2->AzR += fabs(sumZ);

    }
    
    sumX = sumY = sumZ = 0.0;
    for (NI=BegNI(); NI<EndNI(); NI++) {
        if ((degree = NI.GetOutDeg())){
            data1 = GetNDat(NI.GetId());
            
            data1->Ax = DT*((alpha*data1->dxA - beta*data1->dxR)) / (float)degree;
            data1->Ay = DT*((alpha*data1->dyA - beta*data1->dyR)) / (float)degree;
            data1->Az = DT*((alpha*data1->dzA - beta*data1->dzR)) / (float)degree;
            
            data1->AxT = fabs(data1->Ax);
            data1->AyT = fabs(data1->Ay);
            data1->AzT = fabs(data1->Az);
            
            data1->Vx += DT*data1->Ax;
            data1->Vy += DT*data1->Ay;
            data1->Vz += DT*data1->Az;

            data1->x += data1->Vx*DT;
            data1->y += data1->Vy*DT;
            data1->z += data1->Vz*DT;
        }
    }
}

void TParticleNet::RunByStep2(){
    float sumX, sumY, sumZ, r;
    float diffx, diffy, diffz;
    TParticleNet::TNodeI NI, NN;
    TParticleNet::TEdgeI EI;
    PParticle data1, data2;
    int id1, id2, degree;
    
    for (NI=BegNI(); NI<EndNI(); NI++) {
        data1 = GetNDat(NI.GetId());
        data1->dxA = 0;
        data1->dyA = 0;
        data1->dzA = 0;
        data1->dxR = 0;
        data1->dyR = 0;
        data1->dzR = 0;
        data1->AxA = 0;
        data1->AyA = 0;
        data1->AzA = 0;
        data1->AxR = 0;
        data1->AyR = 0;
        data1->AzR = 0;
        
    }
    
    
    // Attraction
    for (EI=BegEI() ; EI < EndEI() ; EI++){
        id1 = EI.GetSrcNId();
        id2 = EI.GetDstNId();
        // Loading node's data;
        data1 = GetNDat(id1);
        data2 = GetNDat(id2);
        diffx = data1->x - data2->x;
        diffy = data1->y - data2->y;
        diffz = data1->z - data2->z;
        r = pow(diffx,2) + pow(diffy,2) + pow(diffz,2);
        r = sqrt(r);
        //        r = 1.0 / (1.0 + exp(-r));
        //        if (r==0) r = 0.00001; // to avoid division per zero when two particles are in the same position
        //        r = 1.0;
        //        sumX = diffx*r;
        //        sumY = diffy*r;
        //        sumZ = diffz*r;
        //        sumX = diffx/r;
        //        sumY = diffy/r;
        //        sumZ = diffz/r;
        sumX = (diffx/r)*(1.0 + 0.0/r);
        sumY = (diffy/r)*(1.0 + 0.0/r);
        sumZ = (diffz/r)*(1.0 + 0.0/r);
        data1->dxA -= sumX;
        data1->dyA -= sumY;
        data1->dzA -= sumZ;
        data2->dxA += sumX;
        data2->dyA += sumY;
        data2->dzA += sumZ;
        
        data1->AxA -= fabs(sumX);
        data1->AyA -= fabs(sumY);
        data1->AzA -= fabs(sumZ);
        data2->AxA += fabs(sumX);
        data2->AyA += fabs(sumY);
        data2->AzA += fabs(sumZ);
    }
    
    // Repulsion
    for (NI=BegNI() ; NI < EndNI() ; NI ++){
        data1 = GetNDat(NI.GetId());
        
        for (NN=NI ; NN < EndNI() ; NN++){
            if (NI.IsNbrNId(NN.GetId()) || NI==NN) continue;
            data2 = GetNDat(NN.GetId());
            diffx = data1->x - data2->x;
            diffy = data1->y - data2->y;
            diffz = data1->z - data2->z;
            r = pow(diffx,2) + pow(diffy,2) + pow(diffz,2);
            //            r = sqrt(r);
            //            if (r==0) r = 0.00001;
            //            sumX = exp(-r)*diffx/(r);
            //            sumY = exp(-r)*diffy/(r);
            //            sumZ = exp(-r)*diffz/(r);
            sumX = diffx/(r);
            sumY = diffy/(r);
            sumZ = diffz/(r);
            data1->dxR -= sumX;
            data1->dyR -= sumY;
            data1->dzR -= sumZ;
            data2->dxR += sumX;
            data2->dyR += sumY;
            data2->dzR += sumZ;
            
            data1->AxR -= fabs(sumX);
            data1->AyR -= fabs(sumY);
            data1->AzR -= fabs(sumZ);
            data2->AxR += fabs(sumX);
            data2->AyR += fabs(sumY);
            data2->AzR += fabs(sumZ);
            
        }
    }
    
    sumX = sumY = sumZ = 0.0;
    for (NI=BegNI(); NI<EndNI(); NI++) {
        if ((degree = NI.GetOutDeg())){
            data1 = GetNDat(NI.GetId());
            
            data1->Ax = DT*((alpha*data1->dxA - beta*data1->dxR)) / (float)degree;
            data1->Ay = DT*((alpha*data1->dyA - beta*data1->dyR)) / (float)degree;
            data1->Az = DT*((alpha*data1->dzA - beta*data1->dzR)) / (float)degree;
            
            data1->AxT = fabs(data1->Ax);
            data1->AyT = fabs(data1->Ay);
            data1->AzT = fabs(data1->Az);
            
            data1->Vx += DT*data1->Ax;
            data1->Vy += DT*data1->Ay;
            data1->Vz += DT*data1->Az;
            
            data1->x += data1->Vx*DT;
            data1->y += data1->Vy*DT;
            data1->z += data1->Vz*DT;
        }
    }
}
*/

void TParticleNet::SaveCentroids(const char *filename){
    ofstream file;

    vector<TCentroid>::iterator i;

    int aux;
    for (aux=1, i=Centroids.begin() ; i!=Centroids.end(); ++i, aux++){
        i->comm_id = aux;
    }

    file.open(filename, ofstream::out);
    file << "% id_centroids; position; # of associated particles; total repulsion; total attraction; sum of dists " << endl;
    for (i = Centroids.begin() ; i != Centroids.end(); ++i){
        for (j=0 ; j<PDIM ; j++) {
            file << i->x[j] << "\t";
        }
        file << (int) i->comm_id << endl;
    }
    //             << i->nparticles << "\t"
    //             << i->totalR << "\t"
    //             << i->totalA << "\t"
    //             << i->error << endl;
    file.close();
}


void TParticleNet::SaveNetworkFromParticle(const char *filename, float epsilon){
    FILE *stream;
    stream = fopen(filename, "w");
    if (!stream) return;
    int i,j,count,N;
    int id1, id2;
    PParticle data1, data2;
    TParticleNet::TNodeI NI, NN;
    int NId;
    float Deg=0.0, value;
    float **adjmat;
    
    N = GetNodes();
    
    adjmat = new float*[N];
    for (i=0 ; i<N ; i++) adjmat[i] = new float[N];
    
    count = 0;
    for (NI=BegNI(); NI<EndNI(); NI++) {
        NId = NI.GetId();
        Deg += GetNI(NId).GetDeg();
        data1 = GetNDat(NI.GetId());
        data1->idx = count++;
    }
    Deg /= (float)(2*N);
    cout << "Original Degree: " << Deg << endl;
    
    
    for (NI=BegNI() ; NI < EndNI() ; NI ++){
        data1 = GetNDat(NI.GetId());
        for (NN=NI ; NN < EndNI() ; NN++){
            if (NI.IsNbrNId(NN.GetId()) || NI==NN) {
                adjmat[data1->idx][data1->idx] = 0;
                continue;
            }
            data2 = GetNDat(NN.GetId());
            value = calcDistance(data1->x,data2->x);
            adjmat[data1->idx][data2->idx] = value;
            adjmat[data2->idx][data1->idx] = value;
        }
    }
    
    Deg = 0.0;
    for (i=0 ; i<N-1 ; i++){
        for (j=i+1 ; j<N ; j++){
            if (adjmat[i][j] < epsilon) {
                fprintf(stream,"%d %d\n",i+1,j+1);
                Deg++;
            }
        }
    }
    Deg /= (float)N;
    Deg *=2.0;
    fclose(stream);
    cout << "Degree: " << Deg << endl;
}

int TParticleNet::Infomap(int &max){
    float Q;
    TCnComV com;
    int i;
    TParticleNet::TEdgeI EI;
    int id1, id2;
  
    PUNGraph G = TUNGraph::New(); 
    
    for (EI=BegEI() ; EI < EndEI() ; EI++){    
        id1 = EI.GetSrcNId();
        id2 = EI.GetDstNId();
        if (!G->IsNode(id1)) { G->AddNode(id1); }
        if (!G->IsNode(id2)) { G->AddNode(id2); }
        G->AddEdge(id1,id2);
    }

    Q = TSnap::Infomap(G, com);

    max = 0;
    for (i=0 ; i<com.Len() ; i++){
        if (max < com[i].Len()) max = com[i].Len();
    }
    cout << "INFOMAP #com: " << com.Len() << " Max: " << max << endl;
 
    return com.Len();
}



/*
 int TParticleNet::RunModel(int maxIT, float minDR, bool verbose){
 TParticleNet::TNodeI NI, NN;
 TParticleNet::TEdgeI EI;
 PParticle data1, data2;
 int steps=0;
 float value;
 float dg1, dg2;
 float tempR;
 float diffX=0.0, dTemp=0.0;
 
 
 if (verbose) cout << " Step#   \\theta_R   \\Delta\\theta_R \n";
 
 for (NI=BegNI(); NI<EndNI(); NI++) {
 data1 = GetNDat(NI.GetId());
 for (i=0 ; i<PDIM ; i++) {
 data1->xO[i] = data1->x[i];
 }
 }
 
 
 do {
 
 for (NI=BegNI(); NI<EndNI(); NI++) {
 data1 = GetNDat(NI.GetId());
 data1->DistF = 0;
 for (i=0 ; i<PDIM ; i++) {
 data1->dA[i] = 0;
 data1->dR[i] = 0;
 data1->xOL[i] = data1->x[i];
 }
 }
 
 // Algorithm CORE (O(n^2))
 
 
 // Attraction O(m))
 for (EI=BegEI() ; EI < EndEI() ; EI++){
 id1 = EI.GetSrcNId();
 id2 = EI.GetDstNId();
 // Loading node's data;
 data1 = GetNDat(id1);
 data2 = GetNDat(id2);
 
 for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - data2->x[i];
 r = 0;
 for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
 r = sqrt(r);
 data1->DistF += r;
 data2->DistF += r;
 if (r<1.0) r=1.0;
 for (i=0 ; i<PDIM ; i++) sum[i] = diff[i]/r;
 for (i=0 ; i<PDIM ; i++) {
 data1->dA[i] -= sum[i];
 data2->dA[i] += sum[i];
 }
 }
 
 // repulsion O(n^2) ... n-k \approx n p/ k small
 for (NI=BegNI() ; NI < EndNI() ; NI ++){
 data1 = GetNDat(NI.GetId());
 for (NN=NI ; NN < EndNI() ; NN++){
 if (NI.IsNbrNId(NN.GetId()) || NI==NN) continue;
 data2 = GetNDat(NN.GetId());
 for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - data2->x[i];
 r = 0;
 for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
 r = sqrt(r);
 if (r<0.001) r = 0.001; // to avoid division per zero when the particles are in the same position
 for (i=0 ; i<PDIM ; i++) sum[i] = exp(-r)*diff[i]/(r);
 
 for (i=0 ; i<PDIM ; i++) {
 data1->dR[i] -= sum[i];
 data2->dR[i] += sum[i];
 }
 }
 }
 
 oldRR2 = oldRR;
 oldRR = RR;
 RR = 0.0;
 for (NI=BegNI(); NI<EndNI(); NI++) {
 if ((degree = NI.GetOutDeg())){
 data1 = GetNDat(NI.GetId());
 tempR = 0;
 for (i=0 ; i<PDIM ; i++) {
 data1->x[i] += DT*((alpha*data1->dA[i] - beta*data1->dR[i]))/degree;
 data1->dAO[i] = ((alpha*data1->dA[i] - beta*data1->dR[i]))/degree;
 tempR += fabs(data1->dA[i]);
 }
 RR += tempR/degree;
 }
 dTemp = 0.0;
 for (i=0 ; i<PDIM ; i++) dTemp += pow(data1->x[i] - data1->xOL[i],2);
 diffX += sqrt(dTemp);
 }
 // end Algorithm CORE
 
 RR = RR*0.1 + oldRR*0.9; // used to make the evolution of RR more stable (due to the coarse numerical integration step DT=1.0).
 ++steps;
 
 value = fabs(oldRR2-RR); // / (float)Particles.size();
 if (verbose) cout << steps << " " << RR << " " << value << endl; //" DiffX " << diffX << endl;
 
 if (value < minDR) { // to avoid oscillations (coarse numerical integration)
 if (steps%2) maxIT = steps;
 else maxIT = steps+1;
 }
 } while (steps<maxIT);
 
 
 diffX = 0.0;
 for (NI=BegNI(); NI<EndNI(); NI++) {
 data1 = GetNDat(NI.GetId());
 dTemp = 0.0;
 for (i=0 ; i<PDIM ; i++) dTemp += pow(data1->x[i] - data1->xO[i],2);
 diffX += sqrt(dTemp);
 }
 DiffX = diffX;
 //    cout << "Diff X: " << diffX << endl;
 return steps;
 }
*/

