// dy/dx:   rk4
// dy/dt:   euler
//Improvcement on 6.8's program

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

//second differentiation
double second_diff_euler(double y[],int n){
    return (y[n+1]-2*y[n]+y[n-1])/(dl*dl);
}
//fourth differentiation
double fourth_diff_euler(double y[],int n){
    if(n==1)return (y[n+2]-4*y[n+1]+5*y[n]-4*y[n-1])    / (dl*dl*dl*dl);
    if(n==N-1)return (4*y[n+1]+5*y[n]-4*y[n-1]+y[n-2])  / (dl*dl*dl*dl);
    return (y[n+2]-4*y[n+1]+6*y[n]-4*y[n-1]+y[n-2])     / (dl*dl*dl*dl);
}


// Calculating the acceleration using y0 array, acceleration store in a0 array
void cal_a_simple(double y0[], double a0[]){
    double ypp=0;
    for(int i=1;i<=N-1;i++){
        ypp=second_diff_euler(y0,i);
        a0[i]=(c*c)*(ypp);
        //printf("%3.8f ", a[i]);
    }
}

// Calculating the acceleration with stiffness using y0 array, acceleration store in a0 array 
void cal_a_with_stiffness(double y0[], double a0[]){
    double ypp=0, ypppp=0;
    for(int i=1;i<=N-1;i++){
        ypp=second_diff_euler(y0,i);
        ypppp=fourth_diff_euler(y0,i);
        a0[i]=(c*c)*(ypp-kappa*ypppp);
    }
    output_file(fp1, a0);

}

void cal_a_with_stf_visco(double y0[], double v0[], double a0[]){
    double ypp=0, ypppp=0, vpp=0, aaa=0;
    for(int i=1;i<=N-1;i++){
        ypp     = second_diff_euler(y0,i);
        ypppp   = fourth_diff_euler(y0,i);
        vpp     = second_diff_euler(v0,i);
        a0[i]=(c*c)* ypp- (kappa*kappa)*ypppp + 2 * b2 * vpp;
        a0[i] *= dt * dt;
        //aaa=a0[i];
    }
    output_file(fp1, a0);

}


// Suppose y_rk4[0][MAXN] is the vertical y, and time t is the horizontal 
// y_rk4[j][MAXN] is the j's derivative of t with respective of time t
// y[MAXN] is the y_value array at time t before take any step
// v[MAXN] is the velocity array at time t before take any step 
// y_rk4[0][MAXN] is the y_value array after stepped
// y_rk4[1][MAXN] is the velocity array after stepped 
double y_rk4[3][MAXN];

// slope_rk4[1][0][MAXN] is slope for velocity( namely, y_rk4[1]) to take the first step
// ...
// slope_rk4[1][3][MAXN] is slope for velocity( namely, y_rk4[1]) to take the fourth step
// slope_rk4[0][0][MAXN] is slope for y( namely, y_rk4[0]) to take the first step
// ...
// slope_rk4[0][3][MAXN] is slope for y( namely, y_rk4[0]) to take the fourth step
double slope_rk4[3][5][MAXN];


// the k th step
// the J th derivative of y: y(j)
// the slope stored in slope_stored
void find_slope_rk4(int j, int k, double slope_stored[]){
    if(j==1){
        // Using y( namely, y_rk4[0]) to calculate acceleration
        if(k==0) cal_a_with_stf_visco(y, v, slope_stored);
        else cal_a_with_stf_visco(y_rk4[0], y_rk4[1], slope_stored);  
    }     
    if(j==0){
        if(k==0){
            for(int i=1;i<N;i++){
                slope_stored[i] = v[i];
            }
        }
        else{
            for(int i=1;i<N;i++){
                slope_stored[i] = y_rk4[1][i];
            }
        }
    } 
}

// J th order ODE
// the j th derivative of y: y(j)
// the k th step's step length=step_len
void take_the_step(int J, int k, double step_len){
    for(int i=1;i<N;i++){
        y_rk4[0][i] = y[i] + slope_rk4[0][k][i] * step_len;
    }
    for(int i=1;i<N;i++){
        y_rk4[1][i] = v[i] + slope_rk4[1][k][i] * step_len;
    }
    
}

// the j th derivative of y: y(j)
// the k th step's step length=step_len
// initial position array before step, final position array after step
void step_rk4(int J, int k, double step_len){
    for(int j=J-1; j>=0;j--){
        find_slope_rk4(j, k, slope_rk4[j][k]);
    }    

    if(k!=3)take_the_step(J, k, step_len);
}

void last_step_rk4(int J, double step_len){
    double delta=0;
    for(int i=0;i<=N;i++){
        delta = slope_rk4[1][0][i] + 2* slope_rk4[1][1][i] + 2 * slope_rk4[1][2][i] + slope_rk4[1][3][i] ;
        delta /= 6;
        v[i] +=  delta * step_len;
    }
    for(int i=0;i<=N;i++){
        delta = slope_rk4[0][0][i] + 2* slope_rk4[0][1][i] + 2 * slope_rk4[0][2][i] + slope_rk4[0][3][i] ;
        delta /= 6;
        y[i] +=  delta * step_len;
    }

}


// J th order ODE
void rk_4(int J,   double step_len){
    step_rk4(2, 0, step_len/2 );
    step_rk4(2, 1, step_len/2 );
    step_rk4(2, 2, step_len/2 );
    step_rk4(2, 3, step_len   ); 
    last_step_rk4(2, step_len);
}

void cal_vx_rk4(){
    while(T<Tmax){
        //output_file(fp1, y);
        rk_4(2, dt);
        T += dt;
    }
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
    while(T<Tmax){
        output_file(fp1, y);
        cal_a_with_stiffness(y,a);
            //the distance cancelled out
        for(int i=1;i<=N-1;i++){
            v[i] += a[i] * dt;
            y[i] += v[i] * dt;
        }
        T += dt;
    }
}


int main(int argc, const char *argv[]){
    fp1=fopen("C:\\Users\\Lenovo\\session_1\\SRS\\23.6\\6.27\\2.csv", "w+");
    
    init();
    //y[N/2]+=0.05;
    init_parabola(0.3);
    cal_vx_rk4();

    fclose(fp1);
    return 0;
    
}