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
    int n=100;
    vector<vector<double>> phi(n+1,vector<double>(n+1,0));
    double del = 2.0/n;
    vector<vector<double>> exact(n+1,vector<double>(n+1,0));
    double start = omp_get_wtime();
    #pragma omp parallel for collapse(2) num_threads(thread_cnt)
    for(int i=0; i<=n; ++i){
        for(int j=0; j<=n; ++j){
            exact[i][j] = f(-1+del*i,-1+del*j);
        }
    }
    double error = 1;
    int i,j,it;
    it = 0;
    while(error>=0.01){
        error = 0;
        #pragma omp parallel for collapse(2) num_threads(thread_cnt) reduction(max:error)
        for(i=1; i<n; ++i){
            for(j=1; j<n; ++j){
                if((i+j)%2==0){
                    phi[i][j] = (phi[i+1][j]+phi[i][j+1]+phi[i-1][j]+phi[i][j-1]+del*del*q(-1+del*i,-1+del*j))/4;
                    error = max(error, abs(phi[i][j]-exact[i][j])/exact[i][j]);
                }
            }
        }
        #pragma omp parallel for collapse(2) num_threads(thread_cnt) reduction(max:error)
        for(i=1; i<n; ++i){
            for(j=1; j<n; ++j){
                if((i+j)%2){
                    phi[i][j] = (phi[i+1][j]+phi[i][j+1]+phi[i-1][j]+phi[i][j-1]+del*del*q(-1+del*i,-1+del*j))/4;
                    error = max(error, abs(phi[i][j]-exact[i][j])/exact[i][j]);
                }
            }
        }
        it++;
        //cout << "Iter no: " << it << ", Error: " << error << "\n";
    }
    double end = omp_get_wtime();
    double unit = omp_get_wtick();
    double time = (end - start);
    cout << "It took " << it << " iterations and " << time << "s to get solution with maximum 1% error of exact solution\n";
    // for(int i=0; i<=n; ++i) cout << phi[i][15] << " ";
    // cout << "\n";
    return 0;
}