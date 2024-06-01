#include<bits/stdc++.h>
#include<omp.h>
using namespace std;

double f(double x, double y){
    return (x*x-1)*(y*y-1);
}

double q(double x, double y){
    return 2*(2-x*x-y*y);
}

int main(int argc, char** argv){
    int thread_cnt;
    if(argc==2) thread_cnt = stoi(argv[1]);
    else{
        cout << "Enter number of threads: ";
        cin >> thread_cnt;
    }
    int n=400;
    vector<vector<double>> phi(n+1,vector<double>(n+1,0));
    double del = 2.0/n;
    vector<vector<double>> exact(n+1,vector<double>(n+1,0));
    double start = omp_get_wtime();
    for(int i=0; i<=n; ++i){
        for(int j=0; j<=n; ++j){
            exact[i][j] = f(-1+del*i,-1+del*j);
        }
    }
    double error = 1;
    int l,i,it;
    it = 0;
    int nd = 2*n - 1;
    while(error>=0.01){
        error = 0;
        int istart,iend;
        for(l=1; l<nd; ++l){
            if(l<n){
                istart = 1;
                iend = l;
            }
            else{
                istart = l-n+2;
                iend = n-1;
            }
            #pragma omp parallel for num_threads(thread_cnt) reduction(max:error)
            for(i=istart; i<=iend; ++i){
                int j = l-i+1;
                phi[i][j] = (phi[i+1][j]+phi[i][j+1]+phi[i-1][j]+phi[i][j-1]+del*del*q(-1+del*i,-1+del*j))/4;
                error = max(error, abs(phi[i][j]-exact[i][j])/exact[i][j]);
            }
        }
        it++;
        cout << "Iter no: " << it << ", Error: " << error << "\n";
    }
    double end = omp_get_wtime();
    double unit = omp_get_wtick();
    double time = (end - start);
    cout << "It took " << it << " iterations and " << time << "s to get solution with maximum 1% error of exact solution\n";
    // for(int i=0; i<=n; ++i) cout << phi[i][15] << " ";
    // cout << "\n";
    return 0;
}