// dy/dx:   rk4
// dy/dt:   euler
//Improvcement on 6.8's program

#include<stdio.h>
#include<math.h>
#define N 127
#define pi 3.14159265
#define MAXN 1024
#define MOD 10000
// N: total number of balls, a,v,x[] the a,v,x of each balls
// originally the balls are uniformly distributed 
// dt: the delta t, T current time, L total length of the spring
double a[MAXN], v[MAXN], x[MAXN], y[MAXN], f[MAXN];
double dt=0.01 / MOD, Tmax=10, T=0, L=1.5, b=0.01, c=100, kappa=0.00001, dl;

FILE *fp1=NULL, *fp2=NULL, *fp3=NULL;
void output_file(FILE *fp0, double array[]);

//二阶差分
double second_diff_euler(double y[],int n){
    return y[n+1]-2*y[n]+y[n-1];
}
//四阶差分
double fourth_diff_euler(double y[],int n){
    if(n==1)return y[n+2]-4*y[n+1]+5*y[n]-4*y[n-1];
    if(n==N-1)return 4*y[n+1]+5*y[n]-4*y[n-1]+y[n-2];
    return y[n+2]-4*y[n+1]+6*y[n]-4*y[n-1]+y[n-2];
}


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

// Calculating the acceleration using y0 array, acceleration store in a0 array
void cal_a_simple(double y0[], double a0[]){
    double ypp=0;
    for(int i=1;i<=N-1;i++){
        ypp=second_diff_euler(y0,i);
        a0[i]=(c*c)/(dl*dl)*(ypp);
        //printf("%3.8f ", a[i]);
    }
}

// Calculating the acceleration with stiffness using y0 array, acceleration store in a0 array 
void cal_a_with_stiffness(double y0[], double a0[]){
    double ypp=0, ypppp=0;
    for(int i=1;i<=N-1;i++){
        ypp=second_diff_euler(y0,i);
        ypppp=fourth_diff_euler(y0,i);
        a0[i]=(c*c)/(dl*dl)*(ypp-kappa*ypppp);
    }
    //output_file(fp1, a0);

}

void DFT_simple(){
    for(int j=0;j<=N;j++){
        f[j]=0;
        for(int i=0;i<=N;i++){
            f[j]+=y[i]*sin(i*pi*j/N);
        }
    }
}

void output(){
    for(int i=0;i<=N;i++){
        printf("%3.8f ", y[i]);
    }
}

void output_file(FILE *fp0, double array[]){
    if(T==0){
        for(int i=0;i<=N;i++)
        fprintf(fp0, "%3.8f,", 0);
        fprintf(fp0, "0\n");
    }
    fprintf(fp0, "%3.8f,", T);
    for(int i=0;i<N;i++)
        fprintf(fp0, "%3.8f,", array[i]);
    fprintf(fp0, "%3.8f\n", array[N]);
}

void cal_vx_euler(){
    int cnt = 0;
    while(T<Tmax){
        if(cnt % MOD == 0 ){
            output_file(fp1, y);
            cnt=0;
        }
        cal_a_with_stiffness(y,a);
            //the distance cancelled out
        for(int i=1;i<=N-1;i++){
            v[i] += a[i] * dt;
            y[i] += v[i] * dt;
        }
        T += dt;
        cnt++;
 
    }
}


int main(int argc, const char *argv[]){
    fp1=fopen("C:\\Users\\Lenovo\\session_1\\SRS\\23.6\\6.27\\datas\\e4.csv", "w+");
    
    init();
    //y[N/2]+=0.05;
    init_parabola(0.3);
    cal_vx_euler();

    fclose(fp1);
    return 0;
    
}