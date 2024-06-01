#include<bits/stdc++.h>
using namespace std;
const double PI = 3.14159265;
double u0(double x){
    if((x>0.5)||(x<0))  return 0;
    else return sin(4*PI*x);
}
int main(){
    double delx = 0.002, delt = 0.0001;
    int n = 2/delx, t = 1/delt; 
    vector<vector<double>> u(t+1, vector<double>(n+1,0)),u_exact(t+1, vector<double>(n+1,0)),u_quick(t+1, vector<double>(n+1,0));
    for(int i=0; i<=n; ++i){
        u[0][i] = u0(i*delx);
        u_exact[0][i] = u0(i*delx);
        u_quick[0][i] = u[0][i];
    }
    for(int i=1; i<=t; ++i){
        for(int j=1; j<n; ++j){
            u[i][j] = u[i-1][j] - (delt/delx)*(u[i-1][j] - u[i-1][j-1]);
            u_exact[i][j] = u0(j*delx - i*delt);
            if (j == 1) u_quick[i][j] = u_quick[i-1][j] - (delt / delx) * (u_quick[i-1][j] - u_quick[i-1][j-1]);
            else u_quick[i][j] = u_quick[i-1][j] - (delt / delx) * ((3.0/8.0) * u_quick[i-1][j] - (7.0/8.0) * u_quick[i-1][j-1] + (1.0/8.0) * u_quick[i-1][j-2] + (3.0/8.0) * u_quick[i-1][j+1]);
        }
    }
    //cout << "Numerical soln at t=0.5:\n ";
    //for(int i=0; i<=n; ++i) cout << u[t/2][i] << " ";
    // cout << "\n\n\nAnalytical soln at t=0.5:\n ";
    // for(int i=0; i<=n; ++i) cout << u_exact[t/2][i] << " ";
    cout << "\n\n\nNumerical soln at t = 1\n";
    for(int i=0; i<=n; ++i) cout << u_quick[t][i] << " ";
    // cout << "\n\n\n Analytical soln at t = 1\n";
    // for(int i=0; i<=n; ++i) cout << u_exact[t][i] << " ";
    //cout << "Numerical soln at t=0 :\n ";
    //for(int i=0; i<=n; ++i) cout << u[0][i] << " ";
    return 0;
}