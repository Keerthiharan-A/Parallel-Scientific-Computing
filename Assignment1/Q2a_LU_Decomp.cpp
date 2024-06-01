#include<bits/stdc++.h>
#include<omp.h>
using namespace std;
double f(double x){
    return sin(5*x);
}
int main(){
    int n = 25;
    double h = 3.0/n;
    vector<double> a(n+1,0), b(n+1,0), c(n+1,0), B(n+1,0), l1(n+1,0), u1(n+1,0);
    B[0] = (-5*f(0)/2.0 + 2*f(h) + f(2*h)/2.0) / h;
    b[0] = 1;
    c[0] = 2;
    u1[0] = 1; 
    for(int i=1; i<n; ++i){
        B[i] = 3.0*(f(h*i+h)-f(h*i-h))/h;
        a[i] = 1;
        b[i] = 4;
        c[i] = 1;
        l1[i] = a[i]/u1[i-1];
        u1[i] = b[i] - l1[i]*c[i-1];
    }
    B[n] = (5*f(n*h)/2.0 - 2*f((n-1)*h) - f((n-2)*h)/2.0) / h;
    b[n] = 1;
    a[n] = 2;
    l1[n] = a[n]/(u1[n-1]);
    u1[n] = b[n] - l1[n]*c[n-1];
    // Solving the LU System
    vector<double> y(n+1,0), x(n+1,0);
    y[0] = B[0];
    for(int i=1; i<=n; ++i)  y[i] = B[i] - l1[i]*y[i-1];
    x[n] = y[n]/u1[n];
    for(int i=n-1; i>=0; --i) x[i] = (y[i] - c[i]*x[i+1])/ u1[i];
    for(int i=0; i<=n; ++i) cout << x[i]  << " ";
    return 0;
}