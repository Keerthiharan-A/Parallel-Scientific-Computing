#include <iostream>
#include <cmath>
#include <vector>
#include <mpi.h>
using namespace std;
const double PI = 3.14159265;

double phi_0(double x) {
    return sin(2 * PI * x);
}
double q(double x, double y) {
    return (x * x + y * y);
}
void Copy(vector<vector<double>>& phi, const vector<vector<double>>& phi_new, int m, int n) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            phi[i][j] = phi_new[i][j];
        }
    }
}
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double start = MPI_Wtime();
    double del = 0.005;
    int n = 2.0 / del;
    int blockSize = n / size;
    vector<vector<double>> phi(blockSize + 2, vector<double>(n + 1, 0)), phi_new(blockSize + 2, vector<double>(n + 1, 0));
    
    // Initialize boundary conditions
    for (int i = 1; i < blockSize + 1; ++i) {
       	phi[i][0] = phi_0(-1 + (rank * blockSize + i - 1) * del);
       	phi_new[i][0] = phi[i][0];
    }
    int iter = 0;
    double error = 1;
    while (true) {
        // Exchange boundary data between neighboring processors
        if (rank > 0) {
            if(rank%2){
            	MPI_Send(&phi[1][0], n + 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            	MPI_Recv(&phi[0][0], n + 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else{
            	MPI_Recv(&phi[0][0], n + 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            	MPI_Send(&phi[1][0], n + 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            }
        }
        if (rank < size - 1) {
            if(rank%2){
            	MPI_Send(&phi[blockSize][0], n + 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            	MPI_Recv(&phi[blockSize + 1][0], n + 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else{
            	MPI_Recv(&phi[blockSize + 1][0], n + 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            	MPI_Send(&phi[blockSize][0], n + 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            }
        }
        error = 0;
        for (int i = 1; i <= blockSize; ++i) {
            if(i == blockSize && rank == size-1) continue;
            double res = 0;
            for (int j = 1; j < n; ++j) {
                phi_new[i][j] = 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1] + del * del * q(-1 + j * del, -1 + (rank * blockSize + i) * del));
                res += abs(phi_new[i][j] - phi[i][j]);
            }
            phi_new[i][n] = (4 * phi_new[i][n-1] - phi_new[i][n-2]) / 3.0;
            res += abs(phi_new[i][n] - phi[i][n]);
            error = max(error,res);
        }
        Copy(phi, phi_new, blockSize+1, n+1);
        iter++;
        //cout << "Iter: " << iter << " Error: " << error << " Rank " << rank << "\n";
        MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        // If error is less than threshold on any processor, break out of the loop
        if (error <= 1e-4) break;
    }
    vector<vector<double>> result(n + 1, vector<double>(n + 1, 0));
    if(rank==0){
    	for (int i = 1; i <= blockSize; ++i) {
    		for (int j = 0; j <= n; ++j) {
    			result[i][j] = phi[i][j];
    		}
    	}
    }
    if (rank != 0) {
    	for (int i = 1; i <= blockSize; ++i) {
        	MPI_Send(&phi[i][0], n + 1, MPI_DOUBLE, 0, rank * blockSize + i, MPI_COMM_WORLD);
    	}
    } else {
    // If rank 0, receive data from other ranks and store it in the 'result' vector
    	for (int src = 1; src < size; ++src) {
    	    for (int i = 1; i <= blockSize; ++i) {
    	        MPI_Recv(&result[src * blockSize + i][0], n + 1, MPI_DOUBLE, src, src * blockSize + i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	    }
    	}
    }
    double end = MPI_Wtime();
    if (rank == 0) {
    	cout << "TIme Taken: " << end - start << "\n"; 
       	cout << "Total iterations: " << iter << endl;
       	cout << "\nPhi at y=0: \n";
       	for (int i = 0; i <= n; ++i) cout << result[n / 2 + 1][i] << " ";
       	cout << "\nPhi at x=0: \n";
       	for (int i = 0; i <= n; ++i) cout << result[i][n / 2] << " ";
       	cout << endl;
    }
    MPI_Finalize();
    return 0;
}
