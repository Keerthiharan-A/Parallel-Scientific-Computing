#include<bits/stdc++.h>
#include<omp.h>
using namespace std;

double f(double x, double y){
    return (x*x-1)*(y*y-1);
}

double q(double x, double y){
    return 2*(2-x*x-y*y);
}

int main(){
    int n=100;
    vector<vector<double>> phi(n+1,vector<double>(n+1,0));
    double del = 2.0/n;
    vector<vector<double>> exact(n+1,vector<double>(n+1,0));
    double start = omp_get_wtime();
    for(int i=0; i<=n; ++i){
        for(int j=0; j<=n; ++j){
            exact[i][j] = f(-1+del*i,-1+del*j);
        }
    }
    double error = -1;
    int it = 0;
    while(abs(error)>=0.01){
        error = 0;
        for(int i=1; i<n; ++i){
            for(int j=1; j<n; ++j){
                phi[i][j] = (phi[i+1][j]+phi[i][j+1]+phi[i-1][j]+phi[i][j-1]+del*del*q(-1+del*i,-1+del*j))/4;
                error = max(error, abs(phi[i][j]-exact[i][j])/exact[i][j]);
            }
        }
        it++;
        //cout << "Iter no: " << it << ", Error: " << error << "\n";
    }
    double end = omp_get_wtime();
    double unit = omp_get_wtick();
    double time = (end - start);
    cout << "It took " << it << " iterations and " << time << "s to get solution with maximum 1% error of exact solution\n";
    // cout << "Gauss Scidel Solution at y=0.5: ";
    // for(int i=0; i<=n; ++i) cout << phi[i][15] << " ";
    // cout << "\n";
    // cout << "Exact Solution at y=0.5: ";
    // for(int i=0; i<=n; ++i) cout << exact[i][15] << " ";
    // cout << "\n";
    return 0;
}