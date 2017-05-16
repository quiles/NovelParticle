

#include "ParticleNet.h"

#include <time.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#define DT 0.01


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
        particle = GetNDat(id);
        particle->indexReal = com;
        if (com>maxCom) maxCom = com;
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
    
//    PParticleNet net = new TParticleNet();
    PParticle particle;
    int id1, id2;
    int lixo;


    while (fscanf(stream, "%d %d", &id1, &id2) == 2){
//        cout << id1 << " - " << id2 << endl;
//        if (!net->IsNode(id1)) {
        if (!IsNode(id1)) {
//            cout << "\n\nDebug\n\n";
            particle = new TParticle();
            for (i=0 ; i<PDIM ; i++) particle->x[i] = (float)(rand()%2000 - 1000) / 10000.0;
            particle->index = NULL;
            particle->indexReal = 0;
            particle->cluster_id = 0;
            particle->refreshed = true;
//            net->AddNode(id1,particle);
            AddNode(id1,particle);
        }
        if (!IsNode(id2)) {
//        if (!net->IsNode(id2)) {
            particle = new TParticle();
            for (i=0 ; i<PDIM ; i++) particle->x[i] = (float)(rand()%2000 - 1000) / 10000.0;
            particle->index = NULL;
            particle->indexReal = 0;
            particle->cluster_id = 0;
            particle->refreshed = true;
            AddNode(id2,particle);
//            net->AddNode(id2,particle);
        }
        if (!IsEdge(id1,id2)) {
//        if (!net->IsEdge(id1,id2)) {
//            net->AddEdge(id1, id2);
//            net->AddEdge(id2, id1);
            AddEdge(id1, id2);
            AddEdge(id2, id1);
        }
    }
    fclose(stream);
}


void TParticleNet::ReloadNetwork(const char *filename){
    FILE *stream;
    stream = fopen(filename, "r+");
    if (!stream) return;

    int id1, id2;
    PParticle particle;
    TParticleNet::TNodeI NI;
    TParticleNet::TEdgeI EI;

    
    for (NI=BegNI(); NI<EndNI(); NI++) {
        particle = GetNDat(NI.GetId());
        particle->refreshed = false;
    }
    
    for (EI=BegEI() ; EI < EndEI() ; EI++){
        id1 = EI.GetSrcNId();
        id2 = EI.GetDstNId();
        DelEdge(id1,id2,false);
        DelEdge(id2,id1,false);
    }
    

    while (fscanf(stream, "%d %d", &id1, &id2) == 2){
        if (!IsNode(id1)) {
            particle = new TParticle();
            for (i=0 ; i<PDIM ; i++) particle->x[i] = (float)(rand()%2000 - 1000) / 10000.0;
            particle->index = NULL;
            particle->indexReal = 0;
            particle->cluster_id = 0;
            particle->refreshed = true;
            AddNode(id1,particle);
        }
        else {
            particle = GetNDat(id1);
            particle->refreshed = true;
        }

        if (!IsNode(id2)) {
            particle = new TParticle();
            for (i=0 ; i<PDIM ; i++) particle->x[i] = (float)(rand()%2000 - 1000) / 10000.0;
            particle->index = NULL;
            particle->indexReal = 0;
            particle->cluster_id = 0;
            particle->refreshed = true;
            AddNode(id2,particle);
        }
        else {
            particle = GetNDat(id2);
            particle->refreshed = true;
        }

        if (!IsEdge(id1,id2)) {
            AddEdge(id1, id2);
            AddEdge(id2, id1);  
        }
    }
    fclose(stream);

    

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


void TParticleNet::NewNode(int node_id){
    PParticle particle;
    if (IsNode(node_id)) return;

    particle = new TParticle();
    for (i=0 ; i<PDIM ; i++) particle->x[i] = (float)(rand()%2000 - 1000) / 10000.0;
    particle->index = NULL;
    particle->indexReal = 0;
    particle->cluster_id = 0;
    AddNode(node_id,particle);
}

void TParticleNet::NewLink(int i, int j){
    if (IsNode(i) && IsNode(j))
       if (!IsEdge(i,j)) {
           AddEdge(i,j);
           AddEdge(j,i);
       }
}

void TParticleNet::DeleteNode(int node_id){
    PParticle particle;
    if (IsNode(node_id)) {
        particle = GetNDat(node_id);
        delete particle;
        DelNode(node_id);
    }
}

void TParticleNet::DeleteLink(int i, int j){
    if (IsNode(i) && IsNode(j)){
       if (IsEdge(i,j)) DelEdge(i,j,false);
    }
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


//void TParticleNet::LoadCommunities(const char *filename, int format){
//    ifstream file;
//    string line;
//    stringstream ss;
//    int node_id;
//    int com;
//    vector<TParticle>::iterator it;
//    
//    numCommunities = 0;
//    
//    file.open(filename, ifstream::in);
//
//    
//    switch (format) {
//        case 1: { // LFR
//            int teste=0;
//            cout << "Loading Loading: " << file.is_open()<< endl;
//            while (file >> node_id && file >> com){
//                for (it = Particles.begin() ; it != Particles.end() && it->node_id != node_id ; ++it);
//                it->indexReal = com;
//                teste++;
//                if (com > numCommunities) numCommunities = com;
//            }
//        }; break;
//        case 2: { // SNAP
//            com = 1;
//            while (getline(file, line)){
//                ss.clear();
//                ss << line;
//                while (ss >> node_id){
//                    for (it = Particles.begin() ; it != Particles.end() && it->node_id != node_id ; ++it) {
//                    }
//                    it->indexReal = com;
//                    if (com > numCommunities) numCommunities = com;
//                }
//                com++;
//            }
//        }; break;
//    }
//    
//    file.close();
//}
//
//// still under development -> it aims to reduce the number of steps when new nodes are added.
//// Some iterations performed only on the new nodes can make the convergence faster...
//
//void TParticleNet::RunForNewNodes(int steps){
//    float sumX, sumY, sumZ, r, dist;
//    vector<TParticle>::iterator i,j;
//    TUNGraph::TNodeI nodeIt;
//    int s;
//    
//    
//    for (s=0 ; s<steps ; s++){
//        for (i = Particles.begin() ; i != Particles.end(); ++i){
//            i->dxA = 0.0;
//            i->dyA = 0.0;
//            i->dzA = 0.0;
//            i->dxR = 0.0;
//            i->dyR = 0.0;
//            i->dzR = 0.0;
//        }
//        
//        for (i = Particles.begin() ; i != Particles.end()-1; ++i){
//            nodeIt = Network->GetNI(i->node_id);
//            for (j = i+1 ; j != Particles.end(); ++j){
//                r = pow(i->x - j->x,2) + pow(i->y - j->y,2) + pow(i->z - j->z,2);
//                r = sqrt(r);
//                if (r==0) r = 0.0001;
//                
//                if (nodeIt.IsNbrNId(j->node_id)) { // i and j are connected
//                    sumX = (j->x - i->x)/r;
//                    sumY = (j->y - i->y)/r;
//                    sumZ = (j->z - i->z)/r;
//                    i->dxA += sumX;
//                    i->dyA += sumY;
//                    i->dzA += sumZ;
//                    j->dxA -= sumX;
//                    j->dyA -= sumY;
//                    j->dzA -= sumZ;
//                }
//                else { // otherwise
//                    dist = exp(-r);
//                    sumX = dist*(j->x - i->x)/r;
//                    sumY = dist*(j->y - i->y)/r;
//                    sumZ = dist*(j->z - i->z)/r;
//                    i->dxR += sumX;
//                    i->dyR += sumY;
//                    i->dzR += sumZ;
//                    j->dxR -= sumX;
//                    j->dyR -= sumY;
//                    j->dzR -= sumZ;
//                }
//            }
//        }
//        
//        for (i = Particles.begin() ; i != Particles.end(); ++i){
//            if (i->degree != 0){
//                i->x += ((alpha*i->dxA - beta*i->dxR)) / (float)i->degree;
//                i->y += ((alpha*i->dyA - beta*i->dyR)) / (float)i->degree;
//                i->z += ((alpha*i->dzA - beta*i->dzR)) / (float)i->degree;
//            }
//        }
//    }
//}
//

void TParticleNet::RunByStep(){
    TParticleNet::TNodeI NI, NN;
    TParticleNet::TEdgeI EI;
    PParticle data1, data2;
    
    for (NI=BegNI(); NI<EndNI(); NI++) {
        data1 = GetNDat(NI.GetId());
        for (i=0 ; i<PDIM ; i++) {
            data1->dA[i] = 0;
            data1->dR[i] = 0;
        }
    }
    
    
    // Attraction
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
        if (r<1.0) r=1.0;
        for (i=0 ; i<PDIM ; i++) sum[i] = diff[i]/r;
        for (i=0 ; i<PDIM ; i++) {
            data1->dA[i] -= sum[i];
            data2->dA[i] += sum[i];
        }
    }

    
    // repulsion
    for (NI=BegNI() ; NI < EndNI() ; NI ++){
        data1 = GetNDat(NI.GetId());
        for (NN=NI ; NN < EndNI() ; NN++){
            if (NI.IsNbrNId(NN.GetId()) || NI==NN) continue;
            data2 = GetNDat(NN.GetId());         
            for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - data2->x[i];
            r = 0;
            for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
            r = sqrt(r);
if (r<0.01) r = 0.01;
            for (i=0 ; i<PDIM ; i++) sum[i] = exp(-r)*diff[i]/(r);

            for (i=0 ; i<PDIM ; i++) {
                data1->dR[i] -= sum[i];
                data2->dR[i] += sum[i];
            }
        }
    }
        
    for (NI=BegNI(); NI<EndNI(); NI++) {
        if ((degree = NI.GetDeg())){
            data1 = GetNDat(NI.GetId());            
            for (i=0 ; i<PDIM ; i++) {
//                data1->dA[i] = (alpha*data1->dA[i])/(float)degree;
//                data1->dR[i] = (beta*data1->dR[i])/(float)degree;
//                data1->x[i] += DT*((data1->dA[i] - data1->dR[i]));
                data1->x[i] += ((alpha*data1->dA[i] - beta*data1->dR[i]))/degree;
            }
        }
    }
}

int TParticleNet::RunModelNumerical(int maxIT, float minDR, bool verbose){
    TParticleNet::TNodeI NI, NN;
    TParticleNet::TEdgeI EI;
    PParticle data1, data2;
    int steps=0;
    float value;
    
    if (verbose) cout << " Step#   \\theta_R   \\Delta\\theta_R \n";
    
    do {
        
        for (NI=BegNI(); NI<EndNI(); NI++) {
            data1 = GetNDat(NI.GetId());
            for (i=0 ; i<PDIM ; i++) {
                data1->dA[i] = 0;
                data1->dR[i] = 0;
            }
        }
        
        // Algorithm CORE (O(n^2))
        
        // Attraction
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
            if (r<1.0) r=1.0;
            for (i=0 ; i<PDIM ; i++) sum[i] = diff[i]/r;
            for (i=0 ; i<PDIM ; i++) {
                data1->dA[i] -= sum[i];
                data2->dA[i] += sum[i];
            }
        }
        
        
        // repulsion
        for (NI=BegNI() ; NI < EndNI() ; NI ++){
            data1 = GetNDat(NI.GetId());
            for (NN=NI ; NN < EndNI() ; NN++){
                if (NI.IsNbrNId(NN.GetId()) || NI==NN) continue;
                data2 = GetNDat(NN.GetId());
                for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - data2->x[i];
                r = 0;
                for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
                r = sqrt(r);
//                if (r<0.01) r = 0.01;
                for (i=0 ; i<PDIM ; i++) sum[i] = exp(-r)*diff[i]/(r);
                
                for (i=0 ; i<PDIM ; i++) {
                    data1->dR[i] -= sum[i];
                    data2->dR[i] += sum[i];
                }
            }
        }
        
        oldRR = RR;
        RR = 0.0;
        float tempR;
        float tempA;
        float AA=0.0;
        for (NI=BegNI(); NI<EndNI(); NI++) {
            if ((degree = NI.GetDeg())){
                data1 = GetNDat(NI.GetId());
                tempR = 0;
                tempA = 0;
                for (i=0 ; i<PDIM ; i++) {
                    data1->x[i] += DT*((alpha*data1->dA[i] - beta*data1->dR[i]))/degree;
//                    tempR += fabs(data1->dR[i]);
//                    tempA += fabs(data1->dA[i]);
                    tempR += pow(data1->dR[i],4);
                    tempA += pow(data1->dR[i],4);
                }
                RR += tempR/degree;
                AA += tempA/degree;
            }
        }
        
        
        // end Algorithm CORE
        
        
//        RR = RR*0.1 + oldRR*0.9; // used to make the evolution of RR more stable (due to the numerical integration step DT=1.0.
        ++steps;
        
        value = fabs(oldRR-RR); // / (float)Particles.size();
        cout << steps << "\t" << AA << "\t" << RR << endl;
//        if (verbose) cout << steps << " " << RR << " " << value << endl;
//        if (verbose) cout << steps << " " << value << endl;
    } while (steps<maxIT && value > minDR);
    return steps;
}


int TParticleNet::RunModel(int maxIT, float minDR, bool verbose){
    TParticleNet::TNodeI NI, NN;
    TParticleNet::TEdgeI EI;
    PParticle data1, data2;
    int steps=0;
    float value;

    if (verbose) cout << " Step#   \\theta_R   \\Delta\\theta_R \n";
    
    do {
        
        for (NI=BegNI(); NI<EndNI(); NI++) {
            data1 = GetNDat(NI.GetId());
            for (i=0 ; i<PDIM ; i++) {
                data1->dA[i] = 0;
                data1->dR[i] = 0;
            }
        }

        // Algorithm CORE (O(n^2))

        // Attraction
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
            if (r<1.0) r=1.0;
            for (i=0 ; i<PDIM ; i++) sum[i] = diff[i]/r;
            for (i=0 ; i<PDIM ; i++) {
                data1->dA[i] -= sum[i];
                data2->dA[i] += sum[i];
            }
        }

    
        // repulsion
        for (NI=BegNI() ; NI < EndNI() ; NI ++){
            data1 = GetNDat(NI.GetId());
            for (NN=NI ; NN < EndNI() ; NN++){
                if (NI.IsNbrNId(NN.GetId()) || NI==NN) continue;
                data2 = GetNDat(NN.GetId());         
                for (i=0 ; i<PDIM ; i++) diff[i] = data1->x[i] - data2->x[i];
                r = 0;
                for (i=0 ; i<PDIM ; i++) r += diff[i]*diff[i];
                r = sqrt(r);
if (r<0.01) r = 0.01;
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
	float tempR;
        for (NI=BegNI(); NI<EndNI(); NI++) {
            if ((degree = NI.GetDeg())){
                data1 = GetNDat(NI.GetId());            
                tempR = 0;
                for (i=0 ; i<PDIM ; i++) {
                    data1->x[i] += ((alpha*data1->dA[i] - beta*data1->dR[i]))/degree;
                    tempR += fabs(data1->dA[i]);
                }
                RR += tempR/degree;
            }
        }


        // end Algorithm CORE
        
        
        RR = RR*0.1 + oldRR*0.9; // used to make the evolution of RR more stable (due to the numerical integration step DT=1.0.
        ++steps;
    
        value = fabs(oldRR2-RR); // / (float)Particles.size();
        if (verbose) cout << steps << " " << RR << " " << value << endl;
    } while (steps<maxIT && value > minDR);
    return steps;
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
    }
}

int TParticleNet::CommunityDetection3(){
    float oldAcc;
    float diffError;
    int steps=0;
    
    toAddCentroid = false;
    toRemoveCentroid = false;
    
    do {
        oldAcc = accError;
        mergeCentroids();
        assignCentroids();
        if (toAddCentroid || toRemoveCentroid){
            if (toRemoveCentroid) removeCentroid();
            if (toAddCentroid) addCentroid();
            assignCentroids();
        }
        computeCentroids();
        accError = accError*0.1 + oldAcc*0.9; // used to make the evolution of the centroid error more stable...
        diffError = fabs(oldAcc - accError); // / (float)Centroids.size();
        ++steps;
    } while (diffError>0.01);
    return steps;
}

int TParticleNet::CommunityDetection(int n_comm){
    float oldAcc;
    float diffError;
    int steps=0;
    TCentroid centroid;
    int c;
    
    Centroids.clear();
    nextComId=1;
    for (c=0 ; c<n_comm ; c++){
        cout << "DEBUG X\n";
        for (i=0 ; i<PDIM ; i++) centroid.x[i] = -1.0 + (float)(rand()%1001) / 500.0;
        centroid.comm_id = nextComId++;
        Centroids.push_back(centroid);
    }
    
    for (c=0 ; c<209 ; c++){
        if (c%5 == 0) {
            mergeCentroids();
            assignCentroids();
//                if (toRemoveCentroid) removeCentroid();
            if (centroidsMerged) addCentroid();
            assignCentroids();
        }
        assignCentroids();
        computeCentroids();
    }
    
    cout << "Debug #com: " << n_comm << " " << Centroids.size() << endl;

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

float TParticleNet::getNormFR(){
    TParticleNet::TNodeI NI;
    PParticle data;
    float DR, Total=0.0;
    
    for (NI=BegNI() ; NI<EndNI(); NI++){
        data = GetNDat(NI.GetId());
        DR = 0.0;
        for (i=0 ; i<PDIM ; i++) {
            DR += data->dR[i]*data->dR[i];
        }
        DR = sqrt(DR);
        Total += DR*DR;
    }
    return sqrt(Total);
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


void TParticleNet::SaveParticlePosition(const char *filename){
    ofstream file;
    TParticleNet::TNodeI NI;
    PParticle data;
    float DR, DA; 

    file.open(filename, ofstream::out);
    file << "% Particle's file -> a snapshot of the particle space\n";
    file << "% Simulation Parameters -> alpha: " << alpha << " beta: " << beta << "Particle dimension: " << PDIM << endl;
    file << "% node_id, ground truth community id, assigned community id, node_id, x1, x2, fa(x1), fa(x2), ..., fr(x1), fr(x2), ...\n";

    for (NI=BegNI() ; NI<EndNI(); NI++){
        data = GetNDat(NI.GetId());
        
        if (data->index) {
            file << fixed << setprecision(2) << NI.GetId() << "\t " 
                                             << data->indexReal << "\t" 
                                             << data->index->comm_id << "\t";
            for (i=0 ; i<PDIM ; i++)  file << fixed << setprecision(2) << data->x[i] << "\t";
            DA = DR = 0.0;
            for (i=0 ; i<PDIM ; i++) {
                DA += data->dA[i]*data->dA[i];
                DR += data->dR[i]*data->dR[i];
            }
            DA = sqrt(DA);
            DR = sqrt(DR);
            file << fixed << setprecision(3) << DA << "\t" << DR << "\t";
            file << endl;
        }
        else {
            file << fixed << setprecision(2) << NI.GetId() << "\t " 
                                             << data->indexReal << "\t" 
                                             << "-1" << "\t";
            for (i=0 ; i<PDIM ; i++)  file << fixed << setprecision(2) << data->x[i] << "\t";
            DA = DR = 0.0;
            for (i=0 ; i<PDIM ; i++) {
                DA += data->dA[i]*data->dA[i];
                DR += data->dR[i]*data->dR[i];
            }
            DA = sqrt(DA);
            DR = sqrt(DR);
            file << fixed << setprecision(3) << DA << "\t" << DR;
            file << endl;
        }
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
//    for (i=Particles.begin() ; i!=Particles.end(); ++i){
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

