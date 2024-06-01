#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <omp.h>
#include <iomanip>

using namespace std;

int neighbour_sum(vector<vector<int>>& l, int i, int j){
    int res = 0;
    if(i>0) res += l[i-1][j];
    if(i<(size(l)-1)) res += l[i+1][j];
    if(j>0) res += l[i][j-1];
    if(j<(size(l)-1)) res += l[i][j+1];
    return res;
}

int main(int argc, char* argv[]) {
    int thread_cnt;
    if (argc == 2) thread_cnt = stoi(argv[1]);
    else {
        cout << "Enter number of threads: ";
        cin >> thread_cnt;
    }

    // Parameters
    double a = 1;
    int N = 1000;
    int N_iter = 2000;
    int ln = 10;
    vector<double> r(6, 0.0);
    r[0] = 1.0;
    vector<double> P(6, 0.0);
    vector<vector<double>> Result(N, vector<double>(2, 0));
    // Initialize random number generator
    srand(static_cast<unsigned int>(time(nullptr)));

    // Initialize rates
    for (int i = 1; i <= 6; ++i) {
        r[i] = 2 * pow(a, i - 2);
    }
    
    omp_set_num_threads(thread_cnt);

    double start_time = omp_get_wtime();
    // Parallelize the outer loop
    #pragma omp parallel for firstprivate(P) schedule(static) shared(r, ln, N, Result, N_iter) default(none) 
    for(int iter = 0; iter < N_iter; ++iter){
        vector<vector<int>> L(ln+2, vector<int>(ln+2, 0));
        double t = 0.0;
    for (int m = 0; m < N; ++m) {
        vector<vector<pair<int, int>>> E(6);
        for (int i = 0; i < ln; ++i) {
            for (int j = 0; j < ln; ++j) {
                if (L[i][j] == 0) {
                    E[0].push_back({i, j});
                } else {
                    int s = neighbour_sum(L, i, j);
                    switch (s) {
                        case 0:
                            E[1].push_back({i, j});
                            break;
                        case 1:
                            E[2].push_back({i, j});
                            break;
                        case 2:
                            E[3].push_back({i, j});
                            break;
                        case 3:
                            E[4].push_back({i, j});
                            break;
                        case 4:
                            E[5].push_back({i, j});
                            break;
                    }
                }
            }
        }
        double R = 0.0;
        vector<int> n(6, 0);
        for (int k = 0; k < 6; ++k) {
            n[k] = E[k].size();
            R += n[k] * r[k];
        }
        int rank = omp_get_thread_num();

        // Store results
        #pragma omp critical
        {
            Result[m][0] += t;
            Result[m][1] += static_cast<double> (n[1] + n[2] + n[3] + n[4]) / (ln * ln);
        }
        
        // Calculate probabilities
        P[0] = n[0] * r[0] / R;
        for (int k = 0; k < 5; ++k) {
            P[k + 1] = P[k] + n[k + 1] * r[k + 1] / R;
        }

        double rd = (double)rand() / RAND_MAX;
        // cout << rd << "\n";
        if (rd <= P[0]) {
            int temp = floor(((double)rand() / RAND_MAX) * max(n[0]-1,0));
            int ei = E[0][temp].first; 
            int ej = E[0][temp].second;
            L[ei][ej] = 1;

            // Apply periodic boundary conditions
            if (ei == 1) {
                L[ln][ej] = 1;
            } else if (ei == ln) {
                L[1][ej] = 1;
            }

            if (ej == 1) {
                L[ei][ln] = 1;
            } else if (ej == ln) {
                L[ei][1] = 1;
            }
        } else {
            for (int p = 0; p < 5; ++p) {
                if (rd > P[p] && rd <= P[p + 1]) {
                    int temp = floor(((double)rand() / RAND_MAX) * max(n[p+1]-1,0));
                    int loc = 0;
                    for (int i = 0; i < p; ++i) loc += n[i];
                    int ei = E[p+1][temp].first;
                    int ej = E[p+1][temp].second;
                    L[ei][ej] = 0;

                    // Apply periodic boundary conditions
                    if (ei == 1) {
                        L[ln][ej] = 0;
                    } else if (ei == ln) {
                        L[1][ej] = 0;
                    }

                    if (ej == 1) {
                        L[ei][ln] = 0;
                    } else if (ej == ln) {
                        L[ei][1] = 0;
                    }
                }
            }
        }
        t += 1 / R;
    }
    }
    for(int i=0; i<N; ++i){
        Result[i][0] /= N_iter;
        Result[i][1] /= N_iter;
    }
    double end_time = omp_get_wtime();
    // Calculate duration
    double duration = end_time - start_time;
    // Output the time taken
    cout << "Time taken: " << duration << " seconds" << endl;
    // Write results to file
    ofstream outFile("results.txt");
    if (!outFile) {
        cerr << "Unable to create file!" << endl;
        return 1;
    }
    outFile << "Time: \n";
    for (int j = 0; j < N; ++j) {
        outFile << fixed << std::setprecision(3) << Result[j][0] << endl;
    }
    outFile << "Surface coverage: \n";
    for(int j = 0; j < N; ++j){
        outFile << fixed << std::setprecision(3) << Result[j][1] << endl;
    }
    outFile.close();
    cout << "Execution completed." << endl;
    return 0;
}