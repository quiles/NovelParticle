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

#define PDIM 3

using namespace std;

class TParticleNet;
class TParticle;
typedef TPt<TParticleNet> PParticleNet;
typedef TParticle* PParticle;

class TCentroid{
public:
    float x[PDIM];
    float error;
    int nparticles; // store the number of associated particles
    int comm_id;
    float totalR; // sum of the repulsion inside the cluster
    float totalA; // sum of the attraction inside the cluster
    TCentroid(){
    };
    ~TCentroid(){
    };
};

class TNode{
public:
  PParticle Val;
public:

  TNode(): Val(NULL){}
  TNode(const PParticle& _Val): Val(_Val){}
  ~TNode();
  operator PParticle() const {return Val;}
  void Save(TSOut& SOut) const {}

};

/*
class TLink{
public:
  bool refreshed;
  float weight;
public:
  TLink(): refreshed(true), weight(1.0){}
  TLink(const bool& _ref, const float _wei): refreshed(_ref), weight(_wei){}
  ~TLink() {};
  operator bool() const {return refreshed;}
  operator float() const {return weight;}
  void Save(TSOut& SOut) const {}
};
*/


class TParticle{
public:
//    int node_id;
    //    TUNGraph::TNodeI node;
    int auxi;
    float x[PDIM];
//    int degree;
    float dA[PDIM];
    float dR[PDIM];
    float Vx, Vy, Vz;
    float Ax, Ay, Az;
    float AxR, AyR, AzR;
    float AxA, AyA, AzA;
    float AxT, AyT, AzT;

    
    //	float erro;
    int cluster_id;
    TCentroid *index; // community id obtained by the algorithm / pointer to the centroid
    int indexReal; // real community id
    vector <int> indexRealH; // real community id
    bool refreshed; // indicates whether a node has been updated
    
public:
    TParticle() {
    };
    ~TParticle(){
    }
//    TParticle(const PParticle &pIn){
//        x = pIn.x; y = pIn.y; z = pIn.z;
//        index = pIn.index; indexReal = pIn.indexReal; cluster_id = pIn.cluster_id;
//    }
    TParticle(TSIn& SIn) {};
    void Save(TSOut& SOut) const {};
    friend class TParticleNet;
};


class TParticleNet : public TNodeNet<TNode> {
private:
    // auxiliar variables
    int i,j,k;
    int id1, id2, degree;
    float sum[PDIM], r;
    float diff[PDIM];

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

    void ResetParticles();
    void computeCentroids();
    void assignCentroids();
    void addCentroid();
    void removeCentroid();
    void mergeCentroids();


public:
    
    TParticleNet();
    ~TParticleNet() {};
    
    static PParticleNet LoadFromFile(const char *filename);
    static PParticleNet LoadClique(int nc, int cs);
    void ReloadNetwork(const char *filename);
    void LoadComFile(const char *filename);
    void LoadComFile(int nc);


    void RunByStep();
    void RunByStep2();
    void RunByStepRadial();
    void RunByStepRadial2();

    void SaveParticlePosition(const char *filename);
    int CommunityDetection3();
    void SetModelParameters(float a, float b, float g){
        alpha=a; beta=b; gamma=g; eta=1.0;
    };
    int getNumCommunities() {return Centroids.size();};

    void NewNode(int node_id);
    void NewLink(int i, int j);
    void DeleteNode(int node_id);
    void DeleteLink(int i, int j);
    float NMI();


    int DEGUB_(int i, int j){
        if (IsEdge(i,j)) return 1;//cout << "Link " << i << "-" << j << endl;
        else return 0;
    };    

    friend class TPt<TParticleNet>;
};

#endif /* PHASESEPARATION_H_ */


