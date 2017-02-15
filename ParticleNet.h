/*
 * PhaseSeparation.h
 *
 *  Created on: Set/10/2015
 *      Author: quiles
 */

#ifndef PHASESEPARATION_H_
#define PHASESEPARATION_H_


#include <iostream>
#include <vector>
#include "Snap.h"

#undef min
#undef max

using namespace std;

class TParticleNet;
class TParticle;
typedef TPt<TParticleNet> PParticleNet;
typedef TParticle* PParticle;

class TCentroid{
public:
    float x, y, z;
    float error;
    int nparticles; // store the number of associated particles
    int comm_id;
    float totalR; // sum of the repulsion inside the cluster
    float totalA; // sum of the attraction inside the cluster
};


class TParticle{
public:
//    int node_id;
    //    TUNGraph::TNodeI node;
    float x, y, z;
//    int degree;
    float dxA, dyA, dzA;
    float dxR, dyR, dzR;
    //	float erro;
    int cluster_id;
    TCentroid *index; // community id obtained by the algorithm / pointer to the centroid
    int indexReal; // real community id
    vector <int> indexRealH; // real community id
    
public:
    TParticle() {};
    TParticle(const TParticle &pIn){
        x = pIn.x; y = pIn.y; z = pIn.z;
        index = pIn.index; indexReal = pIn.indexReal; cluster_id = pIn.cluster_id;
    }
    TParticle(TSIn& SIn) {};
    void Save(TSOut& SOut) const {};
    friend class TParticleNet;
};


class TParticleNet : public TNodeNet<TParticle> {
private:
    float alpha; // attraction strength
    float beta; // repulsion strength
    float gamma; // repulsion decai, constant at 1.0
    float eta; // numerical integration step (delta T), constant at 1.0
    float addCentroidThreshold;
    float mergeCentroidThreshold;
    //    int centroidTransient;
    
    vector <TCentroid> Centroids;

    int nextComId;
    int numCommunities;
    vector<int> numCommunitiesH;
    int numClusters;
    
    
    bool toAddCentroid; // flag to add a new centroid
    bool toRemoveCentroid; // flag to remove a centroid
    int centroid2remove;
    int idCentroidMaxError; // store the id of the centroid with the largest error
    float accError;
    float RR, oldRR, oldRR2;

public:
    
    TParticleNet();
    ~TParticleNet() {};
    
    static PParticleNet LoadFromFile(const char *filename);

    void RunByStep();
    void SaveParticlePosition(const char *filename);

    
    friend class TPt<TParticleNet>;
};

#endif /* PHASESEPARATION_H_ */


