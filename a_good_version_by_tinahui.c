//
//  main.c
//  Initial Position Correction
//
//  Created by LIANG LIAO on 2023/7/3.
//

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

static int N = 255;
double L = 1.5;
double dt = 0.0001;
double tmax = 10;
double t = 0;
double c = 200, e = 10, b = 1000;
int position = 128;
double displacement = 0.005;

//time unit: second
//length unit: meter
//c=f*2L=sqrt(T/(ρA)) (wave speed)
//e=EI/(ρA)=Er^2/(4ρ) (stiffness)
//b: viscoelastisity damping

int output_array(double arr[]){
    for(int n = 0; n <= N; n++){
        printf("%2.10f,", arr[n]);
    }
    printf("\n");
    return 0;
}

int initialize(double a[], double v[], double y[], double x[]){
    for (int n = 0; n<=N; n++){
        v[n] = 0;
        y[n] = 0;
        x[n] = L/N * n;
    }
    return 0;
}

int find_pp(double y[], double pp[]){
    for (int n = 1; n <= N-1; n++){
        pp[n] = y[n+1] - 2*y[n] + y[n-1];
    }
    return 0;
}

int find_pppp(double y[], double pppp[]){
    pppp[1] = y[3] - 4*y[2] + 6*y[1] - 4*y[0];
    pppp[N-1] = - 4*y[N] + 6*y[N-1] - 4*y[N-2] + y[N-3];
    for (int n = 2; n <= N-2; n++){
        pppp[n] = y[n+2] - 4*y[n+1] + 6*y[n] - 4*y[n-1] + y[n-2];
    }
    return 0;
}

int pluck(double y[], int position, double displacement, double a[], double ypp[], double ypppp[]){
    double k1 = displacement/position;
    double k2 = displacement/(N - position);
    for(int n = 1; n <= position; n++){
        y[n] = k1 * n;
    }
    for(int n = position; n <= N; n++){
        y[n] = k2 * (N - n);
    }
    for (int i = 0; i <= 5000; i++){
        find_pp(y, ypp);
        find_pppp(y, ypppp);
        for (int n = 1; n<N; n++){
            a[n] = c * c * ypp[n] - e * ypppp[n];
            if (n == position){
                y[n] = y[n];
            }else{
                y[n] += a[n] * 0.004 * 0.004 / 2;
            }
        }
//        if (i % 100 == 0){
//            output_array(a);
//        }
    }
//    output_array(a);
//    output_array(y);
    return 0;
}

int FFT(double y[], double data[]){
    for (int i = 0; i <= N; i++){
        REAL(data,i) = y[i]; IMAG(data,i) = 0.0;
    }
    gsl_fft_complex_radix2_forward (data, 1, N+1);
    for (int i = 0; i < (N+1)/2; i++){
//        printf ("%d %e %e\n", i, REAL(data,i)/sqrt(N+1), IMAG(data,i)/sqrt(N+1));
        printf("%.10f,", IMAG(data,i)); //data[n] = data[N-n], imaginary part is sin
    }
    printf ("\n");
    return 0;
}

int main(int argc, const char * argv[]){
    double a[N+1], v[N+1], x[N+1], y[N+1];
    double ypp[N+1], ypppp[N+1], vpp[N+1];
    double ytem[N+1], vtem[N+1];
    double ky[4][N+1], kv[4][N+1];
    double data[2*(N+1)];
    int i = 0;
    
    //print parameters
    printf("FFT,dt,c,e,b,position\n");
    printf(" ,%.6f,%.2f,%.2f,%.2f,%d\n", dt, c, e, b, position);
    
    initialize(a, v, y, x);
    pluck(y, position, displacement, a, ypp, ypppp);
    
    for(t = 0; t <= tmax; t = t + dt){
        
        //find k1
        find_pp(y, ypp);
        find_pppp(y, ypppp);
        find_pp(v, vpp);
        for (int n = 1; n<N; n++){
            kv[0][n] = dt * (c*c * ypp[n] - e * ypppp[n] + b * vpp[n]);
            ky[0][n] = dt * v[n] ;
        }
        //take a half step
        for (int n = 1; n<N; n++){
            vtem[n] = v[n] + kv[0][n] / 2;
            ytem[n] = y[n] + ky[0][n] / 2;
        }
        //find k2
        find_pp(ytem, ypp);
        find_pppp(ytem, ypppp);
        find_pp(vtem, vpp);
        for (int n = 1; n<N; n++){
            kv[1][n] = dt * (c*c * ypp[n] - e * ypppp[n] + b * vpp[n]);
            ky[1][n] = dt * vtem[n];
        }
        //take a half step
        for (int n = 1; n<N; n++){
            vtem[n] = v[n] + kv[1][n] / 2;
            ytem[n] = y[n] + ky[1][n] / 2;
        }
        //find k3
        find_pp(ytem, ypp);
        find_pppp(ytem, ypppp);
        find_pp(vtem, vpp);
        for (int n = 1; n<N; n++){
            kv[2][n] = dt * (c*c * ypp[n] - e * ypppp[n] + b * vpp[n]);
            ky[2][n] = dt * vtem[n];
        }
        //take a full step
        for (int n = 1; n<N; n++){
            vtem[n] = v[n] + kv[2][n];
            ytem[n] = y[n] + ky[2][n];
        }
        //find k4
        find_pp(ytem, ypp);
        find_pppp(ytem, ypppp);
        find_pp(vtem, vpp);
        for (int n = 1; n<N; n++){
            kv[3][n] = dt * (c*c * ypp[n] - e * ypppp[n] + b * vpp[n]);
            ky[3][n] = dt * vtem[n];
        }
        //find new v,y
        for (int n = 1; n<N; n++){
            v[n] = v[n] + (kv[0][n] + 2*kv[1][n] + 2*kv[2][n] + kv[3][n])/6;
            y[n] = y[n] + (ky[0][n] + 2*ky[1][n] + 2*ky[2][n] + ky[3][n])/6;
        }
        
        if (i % 500 == 0){
//            FFT(y, data);
            output_array(y);
        }
        i++;
    }
    return 0;
}
