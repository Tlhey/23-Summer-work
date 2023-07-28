#include<stdio.h>
#include<math.h>
#define N 127
#define pi 3.14159265
#define MAXN 1024
// N: total number of balls, a,v,x[] the a,v,x of each balls
// originally the balls are uniformly distributed 
// dt: the delta t, T current time, L total length of the spring
double a[MAXN], v[MAXN], x[MAXN], y[MAXN], f[MAXN];
double dt=0.001, Tmax=10, T=0, L=0.63, c=329.6, kappa=1.25, dl;
double b1 = 1.1, b2 = 2.7e-4;

FILE *fp1=NULL, *fp2=NULL, *fp3=NULL;
void output_file(FILE *fp0, double array[]);

void init(){
    dl=L/N;
    for(int i=1;i<=N+5;i++){
        a[i]=0; v[i]=0; x[i]=dl * i; y[i]=0;
    }
}

// With the maxium height H
void init_sin(double H){
    for(int i=0;i<=N;i++){
        y[i]=sin(pi*i/(double)N)*H;
    }
}

// With the maxium height H
void init_parabola(double H){
    for(int i=0;i<=N;i++){
        y[i]=x[i]*(L-x[i])*4/(double)(L*L)*H;
        //printf("%3.8f ",y[i]);
;    }
}
