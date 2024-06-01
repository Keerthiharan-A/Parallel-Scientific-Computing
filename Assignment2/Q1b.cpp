#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>
using namespace std;
const double PI = 3.14159265;

double u0(double x) {
    if ((x > 0.5) || (x < 0)) return 0;
    else return sin(4 * PI * x);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double delx = 0.002, delt = 0.0001;
    int n = 2 / delx, t = 1 / delt;

    // Calculate local number of grid points
    int n_local = n / size;
    int start_index = rank * n_local;

    // Allocate memory for local data
    
    vector<vector<double>> u(t + 1, vector<double>(n_local+2,0)); // Including virtual points
    // Initialize the initial condition
    for (int i = 0; i <= n_local + 1; ++i) {
        double x = (start_index + i - 1) * delx;
        u[0][i] = u0(x);
    }

    // Perform time stepping
    for (int i = 1; i <= t; ++i) {
        if (rank > 0) {
            if(rank%2==0){
            	MPI_Send(&u[i-1][1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            	MPI_Recv(&u[i-1][0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else{
            	MPI_Recv(&u[i-1][0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            	MPI_Send(&u[i-1][1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            }
        }
        if (rank < size - 1) {
            if(rank%2==0){
            	MPI_Send(&u[i-1][n_local], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            	MPI_Recv(&u[i-1][n_local + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else{
            	MPI_Recv(&u[i-1][n_local + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            	MPI_Send(&u[i-1][n_local], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            }
        }

        // Perform Upwind scheme
        for (int j = 1; j <= n_local; ++j) {
            u[i][j] = u[i-1][j] - (delt / delx) * (u[i-1][j] - u[i-1][j - 1]);
        }
    }

    // Gather all data to rank 0
    vector<double> u_global_1(n+1, 0), u_global_2(n+1, 0);
    MPI_Gather(&u[t/2][1], n_local, MPI_DOUBLE, &u_global_1[0], n_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&u[t][1], n_local, MPI_DOUBLE, &u_global_2[0], n_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Output the result from rank 0
    if (rank == 0) {
        cout << "Numerical solution at t = 0.5 using Upwind scheme:\n\n";
        for (int i = 0; i <= n; ++i) {
            cout << u_global_1[i] << " ";
        }
        cout << "\n";
        cout << "Numerical solution at t = 1 using Upwind scheme:\n\n";
        for (int i = 0; i <= n; ++i) {
            cout << u_global_2[i] << " ";
        }
        cout << "\n";
    }

    MPI_Finalize();
    return 0;
}
