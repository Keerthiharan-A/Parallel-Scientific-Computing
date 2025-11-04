#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include<omp.h>
using namespace std;

vector<double> T0(int N, int s, int B, int k) {
    vector<double> Pi(N+1, 0.0); // Initial distribution
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(-1, 1); // Uniform distribution for step {-1, 0, +1}
    #pragma omp parallel for shared(Pi,k) schedule(dynamic)
    for (int j = 0; j < B; ++j) {
        int p = k; // The particle starts in k
        int c = 1; // The particle contribution is 1
        for (int n = 0; n < s; ++n) {
            int step = dis(gen); // Draw a step
            //cout << step << "\n";
            p += step; // Set a new position index
            // Anti-reflective boundary at p = 1
            if (p < 1) {
                p = 1;
                c = -c; // Change the contribution sign
            }
            // Anti-reflective boundary at p = N
            if (p > N) {
                p = N;
                c = -c; // Change the contribution sign
            }
        }
        // Set the particle count in the final position
        #pragma omp atomic
        Pi[p] += c;
    }
    for (int i = 1; i <= N; ++i) {
        Pi[i] = Pi[i]*N/B;
    }
    return Pi;
}

int main(int argc, char* argv[]) {
    int thread_cnt;
    if (argc == 2) thread_cnt = stoi(argv[1]);
    else {
        cout << "Enter number of threads: ";
        cin >> thread_cnt;
    }
    omp_set_num_threads(thread_cnt);
    int N = 50;
    int B = 200000;
    int k = 13;

    int s_values[] = {200, 400, 800};
    double start_time = omp_get_wtime();

    for (int s : s_values) {
        vector<double> result = T0(N, s, B, k);
        // cout << "Temperature Distribution for s = " << s << ":" << endl;
        // for (int i = 0; i <= N; ++i) {
        //     cout << result[i] << " ";
        //}
        cout << endl;
    }
    double end_time = omp_get_wtime();
    // Calculate duration
    double duration = end_time - start_time;
    // Output the time taken
    cout << "Time taken: " << duration << " seconds" << endl;
    return 0;
}
