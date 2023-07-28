#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 255
#define NMAX 1000 
///gcc test_gsl.c -o test `gsl-config --cflags --libs` -lm
double L=1.5;
double e = 0 , e_max=10, de=1;
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
    }
}

void pluck_relaxation(int k, int position, double displacement){
    double k1 = displacement / position;
    double k2 = displacement / (N - position);
    for (int n=1; n<=position; n++)     {y[k][n] = k1 * n;}
    for (int n=position + 1; n<=N; n++) {y[k][n] = k2 * (N - n);}
    for(int i=0;i<=1000;i++){
        find_a(k);
        for (int n = 1; n < N; n++) {
            y[k][n] += a[k][n] * 0.003 * 0.003 / 2;
        }
        y[k][position] = displacement;
    }
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

void evolve_with_time(int k, double Tmax, double dt){
    double t=0;
    while(t<Tmax){
        rk4_step_by_dt(k, dt);
        t+=dt;
    }
}

void normalize(int k) {
    double mag2 = 0.0;
    for (int n = 0; n <= N; n++) {
        mag2 += pow(y[k][n], 2);
    }
    double mag = sqrt(mag2);
    for (int n = 0; n <= N; n++) {
        y[k][n] = y[k][n]/mag * 0.01;
    }
}

void output_file(double y[]){
    for(int n=0; n<N; n++){
        fprintf(fp1, "%.8f,", y[n]);
    }
    fprintf(fp1, "%.8f\n", y[N]);
}

void output(double y[]){
    for(int n=0; n<N; n++){
        printf( "%.8f,", y[n]);
    }
    printf( "%.8f\n", y[N]);
}


double d_f[2][NMAX];
double cal_vector_delta_mag(int k1, int k2, int time){
    double mag2 = 0.0;
    for(int n=0; n<=N; n++){
        d_f[time][n] = y[k1][n] - y[k2][n];
        mag2 += pow(d_f[time][n], 2);
    }
    double mag = sqrt(mag2);
    return mag;
}

int main() {
    fp1 = fopen("csv/4.csv", "w+");
    while(e<e_max){
        double z0, zn, lam;
        double T=20, dt=0.00001;
        init(1);
        pluck_relaxation(1, 20, 0.05);
        normalize(1);
        init(2);
        pluck_relaxation(2, 21, 0.05);
        normalize(2);
        z0 = cal_vector_delta_mag(1, 2, 0);

        evolve_with_time(2, T, dt);
        evolve_with_time(1, T, dt);
        zn = cal_vector_delta_mag(1, 2, 1);

        lam = 1/T * (log(zn) - log(z0));
        fprintf(fp1,"%.30f,%.30f\n", e, lam);
        output(y[1]);
        printf("e:%f lam:%f process:%f\n", e, lam, e/e_max);
        e += de;
    }
    fclose(fp1);
    printf("Success!\n");
    return 0;
}
