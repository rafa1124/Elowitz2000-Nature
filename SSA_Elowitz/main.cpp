/*
 //  main.cpp
 //  SSAexample
 //
 //  Created by Yifei Li on 10/29/18.
 //  Copyright © 2018 Yifei Li. All rights reserved.
 //
 
 
 20100106 masa
 A Sample of Gillespie Algorithm (Direct Method) for Autocatalytic Reaction Cycle with C
 
 Refer to:
 EGillespie, D.T., Exact stochastic simulation of coupled chemical reactions,
 The journal of physical chemistry, 81(25), 2340--2361, 1977
 http://pubs.acs.org/doi/abs/10.1021/j100540a008
 
 
 chemical reactions
 
 Autocatalytic Reaction Cycle
 
 X0 + X1 -- c0 --> X1 + X1
 X1 + X2 -- c1 --> X2 + X2
 X2 + X0 -- c2 --> X0 + X0
 
 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <iostream>
#include <random>
using namespace std;

#define ENDTIME   1000        // end of time
#define TIMESTEP  0.005        // interval of output
#define OUTFILE   "SSA_Elowitz.dat"    // output file
#define N         12       // number of reaction
#define M         6      // number of chemical species
#define alpha     216
#define alpha0    0.216
#define n         2
#define beta      0.3


int x[M];            // population of chemical species
//double c[N];            // reaction rates
double p[N];            // propencities of reaction
int u1[N];            // data structure for updating x[]
int u2[N][M];            // data structure for updating x[]; all zero by default;

void init(int x[], int u1[N], int u2[N][M]){
    
    // population of chemical species
    x[0] = 0;
    x[1] = 30;
    x[2] = 35;
    x[3] = 0;
    x[4] = 0;
    x[5] = 0;
//    // reaction rates
//    c[0] = 0.1;
//    c[1] = 0.1;
//    c[2] = 0.1;
    
    // data structure for updating the number of chemical species
    // u1[i][j]: i=reaction number, j=chemical species
    // each element is corresponding to the element of u2 array
//    for(int i = 0; i <= N; i++ ){
//        u1[i]=i;
//        }
    
    u2[0][0] = 1;
    u2[1][1] = 1;
    u2[2][2] = 1;
    u2[3][0] = -1;
    u2[4][1] = -1;
    u2[5][2] = -1;
    u2[6][3] = 1;
    u2[7][4] = 1;
    u2[8][5] = 1;
    u2[9][3] = -1;
    u2[10][4] = -1;
    u2[11][5] = -1;
    
    u1[0] = 0;
    u1[1] = 1;
    u1[2] = 2;
    u1[3] = 0;
    u1[4] = 1;
    u1[5] = 2;
    u1[6] = 3;
    u1[7] = 4;
    u1[8] = 5;
    u1[9] = 3;
    u1[10] = 4;
    u1[11] = 5;
   // u2[N]={}
    // u2[i][j]: i=reaction number, j=in(de)crement of chemical species
   
//    u2[0][0] = -1;
//    u2[0][1] =  1;
//    u2[1][0] = -1;
//    u2[1][1] =  1;
//    u2[2][0] = -1;
//    u2[2][1] =  1;
    
}

// function for updating the reaction propencities
void update_p(double p[], int x[]){
    p[0] = alpha/(1+pow(x[5],n))+alpha0; //m1 transcription
    p[1] = alpha/(1+pow(x[3],n))+alpha0; //m2 transcription
    p[2] = alpha/(1+pow(x[4],n))+alpha0; //m3 transcription
    p[3] = x[0];
    p[4] = x[1];
    p[5] = x[2];
    p[6] = beta*x[0];
    p[7] = beta*x[1];
    p[8] = beta*x[2];
    p[9] = beta*x[3];
    p[10] = beta*x[4];
    p[11] = beta*x[5];
    
}


/////////////////////////////////////////////////////////////
//
// functions
//
int select_reaction(double p[], int pn, double sum_propencity, double r){
    int reaction = -1;
    double sp = 0.0;
    int i;
    r = r * sum_propencity;
    for(i=0; i<pn; i++){
        sp += p[i];
        if(r < sp){
            reaction = i;
            break;
        }
    }
    cout<<"reaction is"<<reaction<< endl;
    return reaction;
}

void update_x(int x[], int u1[N], int u2[N][M], int reaction){
    
    x[u1[reaction]] += u2[reaction][u1[reaction]];
    
}

// output format function
void output(FILE *out, double t, int x[], int xn){
    static double output_t = 0.0;
    int i;
    if(output_t <= t){
        fprintf(out, "%f", t);
        for(i=0; i<xn; i++){
            fprintf(out, "\t%d", x[i]);
        }
        fprintf(out, "\n");
        output_t += TIMESTEP;
    }
}



double sum(double a[], int nn){
    int i;
    double s=0.0;
    for(i=0; i<nn; i++)
        s += a[i];
    return(s);
}

int main(void){
    
    // initialization
    double sum_propencity = 0.0;    // sum of propencities
    double tau=0.0;            // step of time
    double t=0.0;            // time
    double r;            // random number
    int reaction;            // reaction number selected
    
    init(x, u1, u2);
    
    //srand(SEED);
    
    FILE *out=fopen(OUTFILE, "w");
    
    // main loop
    while(t < ENDTIME){
        
        
        // output
        output(out, t, x, M);
        
        // update propencity
        update_p(p, x);
        sum_propencity = sum(p, N);
        
//        random_device rd;  //Will be used to obtain a seed for the random number engine
//        mt19937 gen0(rd()); //Standard mersenne_twister_engine seeded with rd()
//        uniform_real_distribution<> dis(0, 1);
//        double tausum = dis(gen0);
        // sample tau
        if(sum_propencity > 0){
            tau = -log((double)rand()/INT_MAX) / sum_propencity;
        }else{
            break;
        }
        
        // select reaction
        // get well distributed r between [0, 1];
        random_device rd;  //Will be used to obtain a seed for the random number engine
        mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        uniform_real_distribution<> dis(0, 1);
        r = dis(gen);
        reaction = select_reaction(p, N, sum_propencity, r);
        cout<< r<< endl;
        // update chemical species
        update_x(x, u1, u2, reaction);
        
        // time
        t += tau;
        
    }
    fclose(out);
    
    
    return(0);
}
