#include<bits/stdc++.h>
#include<omp.h>
using namespace std;
double f(double x){
    return sin(5*x);
}
int main(int argc, char** argv){
    int thread_cnt;
    if(argc==2){
        thread_cnt = stoi(argv[1]);
    }
    else{
        cout << "Enter number of threads: ";
        cin >> thread_cnt;
    }
    int n = 1000;
    double h = 3.0/n;
    int N = log2(n)+1;
    vector<vector<double>> a(N+1,vector<double>(n+1,0)),b(N+1,vector<double>(n+1,0)),c(N+1,vector<double>(n+1,0)),y(N+1,vector<double>(n+1,0));
    vector<vector<double>> alpha(N+1,vector<double>(n+1,0)),beta(N+1,vector<double>(n+1,0));
    y[0][0] = (-5*f(0)/2.0 + 2*f(h) + f(2*h)/2.0) / h;
    b[0][0] = 1;
    c[0][0] = 2;
    for(int i=1; i<n; ++i){
        y[0][i] = 3.0*(f(h*i+h)-f(h*i-h))/h;
        a[0][i] = 1;
        b[0][i] = 4;
        c[0][i] = 1;
    }
    y[0][n] = (5*f(n*h)/2.0 - 2*f((n-1)*h) - f((n-2)*h)/2.0) / h;
    b[0][n] = 1;
    a[0][n] = 2;
    int i,k;
    double start = omp_get_wtime();
    vector<double> x(n+1,0);
    #pragma omp parallel default(none) shared(a,b,c,y,alpha,beta,N,n,x) private(i,k)
    {
    for(k=1; k<=N; ++k){
        #pragma omp for
        for(i=0; i<=n; ++i){
            y[k][i] = y[k-1][i]; b[k][i] = b[k-1][i];
            if(i>=pow(2,k-1)){
                alpha[k][i] = -a[k-1][i]/b[k-1][i-pow(2,k-1)];
                y[k][i] += alpha[k][i]*y[k-1][i-pow(2,k-1)];
                b[k][i] += alpha[k][i]*c[k-1][i-pow(2,k-1)];
            }
            if(i+pow(2,k-1)<=n){
                beta[k][i] = -c[k-1][i]/b[k-1][i+pow(2,k-1)];
                y[k][i] += beta[k][i]*y[k-1][i+pow(2,k-1)];
                b[k][i] += beta[k][i]*a[k-1][i+pow(2,k-1)];
            }
            if(i>=pow(2,k)) a[k][i] = alpha[k][i]*a[k-1][i-pow(2,k-1)];
            if(i<=(n-pow(2,k))) c[k][i] = beta[k][i]*c[k-1][i+pow(2,k-1)];
        }
    }
    for(int i=0; i<=n; ++i) x[i] = y[N][i]/b[N][i];
    }
    double end = omp_get_wtime();
    cout << "Time taken for " << thread_cnt << " threads is "<< end - start << "s\n";
    // for(int i=0; i<=n; ++i)    cout << x[i] << " ";
    // cout << "\n";
    return 0;
}