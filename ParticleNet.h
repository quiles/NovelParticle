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

#define MAXDIM 20

using namespace std;

class TParticleNet;
class TParticle;
class TLink;
class TLinkData;

typedef TPt<TParticleNet> PParticleNet;
typedef TParticle* PParticle;
typedef TLinkData* PLinkData;

class TCentroid{
public:
    float x[MAXDIM];
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

class TParticle {
public:
    int auxi;
    int idx; // used to index an adjcency matrix
    float x[MAXDIM];
    float xO[MAXDIM];
    float xOL[MAXDIM];
    float dA[MAXDIM];
    float dR[MAXDIM];
    float dAO[MAXDIM];
    float dRO[MAXDIM];
//    float Vx, Vy, Vz;
//    float Ax, Ay, Az;
//    float AxR, AyR, AzR;
//    float AxA, AyA, AzA;
//    float AxT, AyT, AzT;
    float DistF, DistI;
    float ComCentrality; 
// DBSCAN
    bool visited;
    int neighbours;
    int cluster_db;
    int degree;

    
    //	float erro;
    int cluster_id;
    TCentroid *index; // community id obtained by the algorithm / pointer to the centroid
    int indexReal; // real community id
    vector <int> indexRealH; // real community id
    bool refreshed; // indicates whether a node has been updated
    
public:
    TParticle(int dim) {
	int i;
        for (i=0 ; i<dim ; i++) {
            x[i] = (float)(rand()%2000 - 1000) / 10000.0;
            dA[i] = dR[i] = dAO[i] = dRO[i] = 0.0;
        }
        index = NULL;
        indexReal = 0;
        cluster_id = 0;
        refreshed = true;
    };
    ~TParticle(){
    }
    friend class TParticleNet;
};


class TLink : public TVoid{
public:
  PLinkData Val;
public:

   TLink(): Val(NULL){}
   TLink(const PLinkData _Val): Val(_Val){};
   ~TLink() {};
   operator PLinkData() const {return Val;}
//  void Save(TSOut& SOut) const {}
   friend class TParticleNet;
};


class TLinkData{
public:
    float weight;
    float distance;
    int age;
    bool refreshed;
    TLinkData(): refreshed(true), weight(1.0), age(0), distance(0.0){}
    TLinkData(const bool& _ref, const float _wei, const int _age, const float _dist): refreshed(_ref), weight(_wei), age(_age), distance(_dist){}
    ~TLinkData() {};
    friend class TParticleNet;
};




//class TParticleNet : public TNodeNet<TNode> {
class TParticleNet : public TNodeEDatNet<TNode,TLink> {
private:
    // auxiliar variables
    int i,j,k;
    int id1, id2;
    float netDegree;
    float degree;
    float sum[MAXDIM], r;
    float diff[MAXDIM];
    int PDIM;
    int tag;

    float alpha; // attraction strength
    float beta; // repulsion strength
    float gamma; // repulsion decai, constant at 1.0
    float eta; // numerical integration step (delta T), constant at 1.0
    float addCentroidThreshold;
    float mergeCentroidThreshold;
    bool centroidsMerged;
    //    int centroidTransient;
    float DiffX[4], DiffR[4];
    float VarMass[4];


    vector <TCentroid> Centroids;

    int nextComId;
    int nextComDBId;
    int numCommunities;
    vector<int> numCommunitiesH;
    int numClusters;
    
    
    bool toAddCentroid; // flag to add a new centroid
    bool toRemoveCentroid; // flag to remove a centroid
    int centroid2remove;
    int idCentroidMaxError; // store the id of the centroid with the largest error
    int idCentroidMinError; // store the id of the centroid with the largest error
    float accError;
    float RR, oldRR, oldRR2;

    void ResetParticles();
    void computeCentroids();
    void assignCentroids();
    void addCentroid();
    void resetCentroid();
    void removeCentroid();
    void mergeCentroids();
    
    float calcDistance(float *d1, float *d2){
        int i;
        float res=0;
        for (i=0 ; i<PDIM ; i++) res += (d1[i]-d2[i])*(d1[i]-d2[i]);
        return sqrt(res);
    };


public:
    
    TParticleNet(int dim);
    ~TParticleNet() {};
    
    void LoadFromFile(const char *filename);
    void ReloadNetwork(const char *filename);
    void LoadComFile(const char *filename);

    void RunByStep();
    void RunByStep2();
    void RunByStepRadial();
    void RunByStepRadial2();
    int RunModel(int maxIT, float minDR, bool verbose);
    int RunModel2(int maxIT, float minDR, bool verbose);

    int RunModelNumerical(int maxIT, float minDR, bool verbose);

    void SaveParticlePosition(const char *filename);
    void SaveParticleForces(const char *filename);
    void SaveCentroids(const char *filename);
    void SaveMeasures(const char *filename);
    void SaveNetworkFromParticle(const char *filename, float epsilon);

    int CommunityDetection3();
    void CommunityDetectionDB(float epsilon);
    
    void ChangeNetwork(int op, float per);


    float getDiffX(int d) {return DiffX[d];};
    float getDiffR(int d) {return DiffR[d];};
    float getVar(int d) {return VarMass[d];};
    float getNetDegree() {return netDegree;};
    

    void SetModelParameters(float a, float b, float g){
        alpha=a; beta=b; gamma=g; eta=1.0;
/*
        TParticleNet::TEdgeI EI;
        TLink data;
        int id1, id2;
        for (EI=BegEI() ; EI < EndEI() ; EI++){
            id1 = EI.GetSrcNId();
            id2 = EI.GetDstNId();
            data = GetEDat(id1,id2);
            cout << id1 << "-" << id2 << "[AGE: " << data.age << "] [Weight: " << data.weight << "]\n";
        }
*/
    };
    int getNumCommunities() {return Centroids.size();};
    float getNormFR(int d);
    float printCentroidsError();
    int sizeLargeCom();
    int Infomap(int &max);

    void NewNode(int node_id);
    void NewLink(int i, int j);
    void DeleteNode(int node_id);
    void DeleteLink(int i, int j);
    float NMI();
    void SetTag(int t);


    int DEGUB_(int i, int j){
        if (IsEdge(i,j)) return 1;//cout << "Link " << i << "-" << j << endl;
        else return 0;
    };    

    friend class TPt<TParticleNet>;
};

#endif /* PHASESEPARATION_H_ */


