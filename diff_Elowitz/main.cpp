//
//  main.cpp
//  diff_Elowitz
//
//  Created by Yifei Li on 10/29/18.
//  Copyright © 2018 Yifei Li. All rights reserved.
//

//
//  main.cpp
//  test
//
//  Created by Yifei Li on 10/29/18.
//  Copyright © 2018 Yifei Li. All rights reserved.
// Euler Method sovling differentail equation
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#define OUTFILE2   "diff_Elowitz_tmax500.dat"
double alpha = 216;
double alpha0 = 0.216;
double n =2;
double beta = 0.3;
using namespace std;
/* The initial position of the vibration is 4.
 The initial velocity of the vibration is 0 */

// differential equations group;

double dm1(double m1, double p3)
{
    
    return -m1+alpha/(1+pow(p3,n))+alpha0;
}

double dm2(double m2, double p1)
{
    return -m2+alpha/(1+pow(p1,n))+alpha0;
}

double dm3(double m3, double p2)
{
    return -m3+alpha/(1+pow(p2,n))+alpha0;
}

double dp1(double p1, double m1)
{
    return -beta*(p1-m1);
}

double dp2(double p2, double m2)
{
    return -beta*(p2-m2);
}

double dp3(double p3, double m3)
{
    return -beta*(p3-m3);
}

double nexty(double k, double previousy, double h){
    return previousy+h*k;
}

vector<double> EulerMethodODE(double(*pFun)(double r, double y), double t0, double tmax, double initial_y0, int nSteps, double r){
    vector<double> y; // initialization of the output vector
    y.push_back(initial_y0);
    // double t = t0;
    double h = (tmax-t0)/nSteps;
    for (int i = 1; i <= nSteps; i++) {
        double k = pFun(y[i-1],r);
        y.push_back(y[i-1]+h*k);
    }
    return y;
}

// output format function
void output(FILE *out, int j, double p1, double p2, double p3){
    static double output_i = 0.0;
    if(output_i <= j){
        fprintf(out, "\t%d", j);
        fprintf(out, "\t%f", p1);
        fprintf(out, "\t%f", p2);
        fprintf(out, "\t%f", p3);
        fprintf(out, "\n");
        output_i ++;
    }
}

int main() {
    //    double (*pFun)(double,double) = &dm1;
    //    vector<double> BT = EulerMethodODE(pFun,0,5,100,50,0.05);
    
    // parameters initialization time, Steps;
    double tmax = 500;
    double t0 = 0;
    double nSteps = 5000000;
    // calculate step length;
    double h = (tmax-t0)/nSteps;
    // initial value of differential equation;
    
    double initial_m1=0;
    double initial_m2=35;
    double initial_m3=35;
    double initial_p1=0;
    double initial_p2=0;
    double initial_p3=0;
    
    // initialize result vector;
    FILE *out=fopen(OUTFILE2, "w");
    
    vector<double> m1;
    vector<double> m2;
    vector<double> m3;
    vector<double> p1;
    vector<double> p2;
    vector<double> p3;
    
    m1.push_back(initial_m1);
    m2.push_back(initial_m2);
    m3.push_back(initial_m3);
    p1.push_back(initial_p1);
    p2.push_back(initial_p2);
    p3.push_back(initial_p3);
    
    for (int i = 1; i <= nSteps; i++) {
        
        
        double km1 = dm1(m1[i-1],p3[i-1]);
        double km2 = dm2(m2[i-1],p1[i-1]);
        double km3 = dm3(m3[i-1],p2[i-1]);
        double kp1 = dp1(p1[i-1],m1[i-1]);
        double kp2 = dp2(p2[i-1],m2[i-1]);
        double kp3 = dp3(p3[i-1],m3[i-1]);
        
        output(out,i,nexty(kp1,p1[i-1],h),nexty(kp2,p2[i-1],h),nexty(kp3,p3[i-1],h));
        
        m1.push_back(nexty(km1,m1[i-1],h));
        m2.push_back(nexty(km2,m2[i-1],h));
        m3.push_back(nexty(km3,m3[i-1],h));
        p1.push_back(nexty(kp1,p1[i-1],h));
        p2.push_back(nexty(kp2,p2[i-1],h));
        p3.push_back(nexty(kp3,p3[i-1],h));
    }
    
    fclose(out);
    //    ofstream outfile2("test.txt", ios::out | ofstream::binary);
    //    copy(y.begin(), y.end(), ostream_iterator<double>(outfile2));
    //    outfile2.close();
    cout<<dm1(100,1)<<p1[4900000]<<"<<P2>>"<<p3[4900000]<<"<<P2>>"<<p2[4900000];
    
}
