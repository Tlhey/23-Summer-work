// This is used to test whether the energy conserved 
// and Choosing the optimized dt, N, e, b, 


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 1023
#define NMAX 1100 
double L=1.5;
double e = 0 , e_max=20, de=1;
double c = 250;
double b = 0;
double a[4][NMAX], v[4][NMAX], y[4][NMAX], x[NMAX];
FILE *fp1 = NULL;

void init(int k){
    for(int n=0; n<=N+1; n++){
        a[k][n]=0; v[k][n]=0; y[k][n]=0; x[n]=L/N*n;
    }
}

double ypp, ypppp, vpp;
void find_a(int k) {
    for(int n=1; n<N; n++){
        vpp = v[k][n+1] - 2*v[k][n] + v[k][n-1];
        ypp = y[k][n+1] - 2*y[k][n] + y[k][n-1];
        if(n==1)        ypppp = y[k][n+2]- 4*y[k][n+1]+ 6*y[k][n] - 4*y[k][n-1] ; 
        else if(n==N-1) ypppp =          - 4*y[k][n+1]+ 6*y[k][n] - 4*y[k][n-1] + y[k][n-2];
        else            ypppp = y[k][n+2]- 4*y[k][n+1]+ 6*y[k][n] - 4*y[k][n-1] + y[k][n-2];
        a[k][n] = c*c*ypp - e*ypppp + b*vpp;
        //printf("%d %.3f %.3f %.3f\n", n, ypp, ypppp, vpp);
    }
}

void pluck_relaxation(int k, int position, double displacement){
    //printf("%.5f", displacement);
    double k1 = displacement / position;
    double k2 = displacement / (N - position);
    for (int n=1; n<=position; n++)     {y[k][n] = k1 * n;}
    for (int n=position + 1; n<=N; n++) {y[k][n] = k2 * (N - n);}
    // for(int i=0;i<=10000000;i++){
    //     find_a(k);
    //     for (int n = 1; n < N; n++) {
    //         y[k][n] += a[k][n] * 0.000001;
    //     }
    //     y[k][position] = displacement;
    // }
}

double kv[4][10000];
double ky[4][10000];
void rk4_step_by_dt(int k, double dt){
    find_a(k);
    for (int n = 1; n<N; n++){
        kv[0][n] =  a[k][n];
        ky[0][n] =  v[k][n] ;
    }
    for (int n = 1; n<N; n++){
        v[0][n] = v[k][n] + dt*kv[0][n]/2;
        y[0][n] = y[k][n] + dt*ky[0][n]/2;
    }

    find_a(0);
    for (int n = 1; n<N; n++){
        kv[1][n] =  a[0][n];
        ky[1][n] =  v[0][n] ;
    }
    for (int n = 1; n<N; n++){
        v[0][n] = v[k][0] + dt*kv[1][n]/2;
        y[0][n] = y[k][0] + dt*ky[1][n]/2;
    }
        
    find_a(0);
    for (int n = 1; n<N; n++){
        kv[2][n] =  a[0][n];
        ky[2][n] =  v[0][n] ;
    }
    for (int n = 1; n<N; n++){
        v[0][n] = v[k][0] + dt*kv[2][n];
        y[0][n] = y[k][0] + dt*ky[2][n];
    }
    find_a(0);
    for (int n = 1; n<N; n++){
        kv[3][n] =  a[0][n];
        ky[3][n] =  v[0][n] ;
    }
    for (int n = 1; n<N; n++){
        v[k][n] = v[k][n] + (kv[0][n] + 2*kv[1][n] + 2*kv[2][n] + kv[3][n])*dt/6;
        y[k][n] = y[k][n] + (ky[0][n] + 2*ky[1][n] + 2*ky[2][n] + ky[3][n])*dt/6;
    }
}



void output_file(double y[]){
    for(int n=0; n<N; n++){
        fprintf(fp1, "%.10f,", y[n]);
    }
    fprintf(fp1, "%.10f\n", y[N]);
}

void evolve_with_time(int k, double Tmax, double dt){
    double t=0, energy=0;
    int cnt1=0, cnt2=0;;
    while(t<Tmax){
        energy=0;
        rk4_step_by_dt(k, dt);
        //output_file(a[1]);
        for(int n=0;n<=N;n++){
            energy += v[k][n]*v[k][n];
        }
        if(cnt2==1000){fprintf(fp1, "%.10f,%.30f\n", t, energy);cnt2=0;}
        if(cnt1==1000000){printf("Process:%.8f energy:%.10f\n", t/Tmax, energy);cnt1=0;}
        t+=dt;
        cnt1++; cnt2++;
    }
}

void normalize(int k) {
    double mag2 = 0.0;
    for (int n = 0; n <= N; n++) {
        mag2 += pow(y[k][n], 2);
    }
    double mag = sqrt(mag2);
    for (int n = 0; n <= N; n++) {
        y[k][n] = y[k][n] / mag * 0.01;
    }
}

int main() {
    fp1 = fopen("csv/e2.csv", "w+");

    double T=10000, dt=0.001;
    init(1);
    pluck_relaxation(1, 20, 0.005);
    evolve_with_time(1, T, dt);

    fclose(fp1);
    printf("Success!\n");
    return 0;
}
