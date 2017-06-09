
//#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <iostream>
#include <fstream>


using namespace std;

void MaxSTree(float ***matF, bool ***mat, int N);


void NetworkFromMatrix_E(const char *input, const char *output, float epsilon, bool MST){
    FILE *stIn, *stOut;
    int N,i,j;
    bool **mat;
    float **matF;
    float value;
    float degree=0.0;
    
    if ((stIn = fopen(input,"r")) && (stOut = fopen(output,"w")) ) {
        fscanf(stIn,"%d", &N);
        mat = new bool*[N];
        matF = new float*[N];
        for (i=0 ; i<N ; i++){
            mat[i] = new bool[i+1];
            matF[i] = new float[i+1];
            for (j=0 ; j<i+1 ; j++) mat[i][j] = false;
        }
        for (i=0 ; i<N-1 ; i++){
            for (j=i+1 ; j<N ; j++){
                fscanf(stIn,"%f", &value);
                matF[j][i] = value;
            }
        }
        fclose(stIn);
        
        if (MST) MaxSTree(&matF, &mat, N);
        
        for (i=0 ; i<N-1 ; i++){
            for (j=i+1 ; j<N ; j++){
                if (matF[j][i] >= epsilon) mat[j][i] = true;
            }
        }
        
        for (i=0 ; i<N-1 ; i++){
            for (j=i+1 ; j<N ; j++){
                if (mat[j][i]) {
                    fprintf(stOut,"%d %d\n",i+1,j+1);
                    degree++;
                }
            }
        }
        fclose(stOut);
        
        printf("Degree: %.2f\n", degree*2.0 / (float)N);
        
        for (i=0 ; i<N ; i++) {
            delete[] matF[i];
            delete[] mat[i];
        }
        delete[] matF;
        delete[] mat;
    }
    
}

void NetworkFromMatrix_K(const char *input, const char *output, int K, bool MST){
    FILE *stIn, *stOut;
    int N,i,j,k;
    bool **mat;
    float **matF;
    float value;
    float degree=0.0;

    
    if ((stIn = fopen(input,"r")) && (stOut = fopen(output,"w")) ) {
        fscanf(stIn,"%d", &N);
        mat = new bool*[N];
        matF = new float*[N];
        for (i=0 ; i<N ; i++){
            mat[i] = new bool[i+1];
            matF[i] = new float[i+1];
            for (j=0 ; j<i+1 ; j++) mat[i][j] = false;
        }
        for (i=0 ; i<N-1 ; i++){
            for (j=i+1 ; j<N ; j++){
                fscanf(stIn,"%f", &value);
                matF[j][i] = value;
            }
        }
        fclose(stIn);
        
        if (MST) MaxSTree(&matF, &mat, N);

        float max;
        int pmax;
        
        for (i=0 ; i<N-1 ; i++){
            for (k=0 ; k<K ; k++){
                pmax = i+1;
                max = matF[i+1][i];
                for (j=i+2 ; j<N ; j++){
                    if (matF[j][i] > max){
                        pmax = j;
                        max = matF[j][i];
                    }
                }
                mat[pmax][i] = true;
                matF[pmax][i] = 0.0;
            }
        }
        
        for (i=0 ; i<N-1 ; i++){
            for (j=i+1 ; j<N ; j++){
                if (mat[j][i]) {
                    fprintf(stOut,"%d %d\n",i+1,j+1);
                    degree++;
                }
            }
        }
        fclose(stOut);

        printf("Degree: %.2f\n", degree*2.0 / (float)N);

        
        for (i=0 ; i<N ; i++) {
            delete[] matF[i];
            delete[] mat[i];
        }
        delete[] matF;
        delete[] mat;
    }
}

void MaxSTree(float ***matF, bool ***mat, int N){
    int *nodes, nnodes;
    int i,pi,pj,j, pmax, tmp;
    float max;
    bool *connect;
    
    nnodes = N;
    connect = new bool[N];
    nodes = new int[N];
    for (i=0 ; i<N ; i++) {
        nodes[i] = i;
        connect[i] = false;
    }


    for (i=0 ; i<N-1 ; i++){
        pi = nodes[i];
        pmax = i+1;
        pj = nodes[i+1];
        if (pi>pj) max = (*matF)[pi][pj];
        else max = (*matF)[pj][pi];
//        printf("[%d] -> ",i+1);

        for (j=i+2 ; j<N ; j++){
            pj = nodes[j];

            if (pj>pi){
//                printf("%d(%f) ",j+1,matF[pj][pi]);
                if (max < (*matF)[pj][pi]){
                    max = (*matF)[pj][pi];
                    pmax = j;
                }
            } else {
//                printf("%d(%f) ",j+1,matF[pi][pj]);
                if (max < (*matF)[pi][pj]){
                    max = (*matF)[pi][pj];
                    pmax = j;
                }
            }
        }

//        printf("[%d]\n",pmax);
        if (max > 0) {
            connect[i] = true;
            tmp = nodes[i+1];
            nodes[i+1] = nodes[pmax];
            nodes[pmax] = tmp;
        }
    }
    for (i=0 ; i<N-1 ; i++){
        if (!connect[i]) continue;
        pi = nodes[i];
        pj = nodes[i+1];
        if (pj>pi) (*mat)[pj][pi] = true;
        else (*mat)[pi][pj] = true;
//        printf("%d -> %d\n",pi+1,pj+1);
    }
    

    delete[] nodes;
    delete[] connect;
}

void GenerateNetworks(const char *source, int opt, float ke, bool mst){
    string saveName;
    int it, st, count, itc, i;
    char out[256];
    bool firstIt=true;
    int cInfo, maxInfo;
    string fName, mst_c;
    FILE *stream;
    ifstream file(source);
    count = 1;
    
    if (mst) mst_c = "st";
    else mst_c = "";

    while  (file >> fName){
        cout << fName << endl;
        if (opt == 1) sprintf(out,"ep_%.3f%s.txt",ke,mst_c.c_str());
        else sprintf(out,"K_%d%s.txt",(int)ke,mst_c.c_str());
        saveName = fName;
        saveName.replace(fName.size()-3,3,out);
        
        if (opt==1) NetworkFromMatrix_E(fName.c_str(), saveName.c_str(), ke, mst);
        else NetworkFromMatrix_K(fName.c_str(), saveName.c_str(), (int)ke, mst);
    }
}




