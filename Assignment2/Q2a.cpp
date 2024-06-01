#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
const double PI = 3.14159265;

double phi_0(double x) {
    return sin(2 * PI * x);
}

double q(double x, double y) {
    return (x * x + y * y);
}

double norm(const vector<vector<double>>& a) {
    int n = a.size();
    double res = 0;
    for (int i = 0; i < n; ++i) {
        double temp = 0;
        for (int j = 0; j < n; ++j) temp += a[i][j];
        res = max(res, temp);
    }
    return res;
}

void Copy(vector<vector<double>>& phi, const vector<vector<double>>& phi_new) {
    int n = phi.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            phi[i][j] = phi_new[i][j];
        }
    }
}

int main() {
    double del = 0.1;
    int n = 2.0 / del;
    vector<vector<double>> phi(n + 1, vector<double>(n + 1, 0)), phi_new(n + 1, vector<double>(n + 1, 0)), res(n + 1, vector<double>(n + 1, 0));
    // Initialize boundary conditions
    for (int i = 0; i <= n; ++i) {
        phi[0][i] = phi_0(-1 + i * del);
        phi_new[0][i] = phi[0][i];
    }
    int iter = 0;
    double error = 1;

    while (error > 1e-4) {
        for (int j = 1; j < n; ++j) {
            for (int i = 1; i < n; ++i) {
                phi_new[i][j] = 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1] + del * del * q(-1 + i * del, -1 + j * del));
                res[i][j] = abs(phi_new[i][j] - phi[i][j]);
            }
            phi_new[n][j] = (4 * phi_new[n - 1][j] - phi_new[n - 2][j]) / 3.0;
            res[n][j] = abs(phi_new[n][j] - phi[n][j]);
        }
        error = norm(res);
        cout << iter << " " << error << "\n";
        Copy(phi, phi_new);
        iter++;
    }
    cout << "It took " << iter << " iterations to converge to the error " << error << endl;

    cout << "\n\nPhi at y=0: \n\n";
    for (int i = 0; i <= n; ++i) cout << phi[i][n / 2] << " ";
    cout << "\n\nPhi at x=0: \n\n";
    for (int i = 0; i <= n; ++i) cout << phi[n / 2][i] << " ";
    return 0;
}