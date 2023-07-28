#include<stdio.h>
#include<math.h>
#define N 255
double dt = 0.0001,tmax = 10, t = 0;
double L=1.5, dl;
double c = 200, e = 10, b = 1000;
double displacement;
double a[N+1], v[N+1], x[N+1], y[N+1];
double ypp[N+1], ypppp[N+1], vpp[N+1];
double ytem[N+1], vtem[N+1];

double second_diff_euler(double y[],int n){
    return (y[n+1]-2*y[n]+y[n-1])/(dl*dl);
}
//fourth differentiation
double fourth_diff_euler(double y[],int n){
    if(n==1)return (y[n+2]-4*y[n+1]+5*y[n]-4*y[n-1])    / (dl*dl*dl*dl);
    if(n==N-1)return (4*y[n+1]+5*y[n]-4*y[n-1]+y[n-2])  / (dl*dl*dl*dl);
    return (y[n+2]-4*y[n+1]+6*y[n]-4*y[n-1]+y[n-2])     / (dl*dl*dl*dl);
}

void init(){
    dl = L / N;
}

int main(){
    
}