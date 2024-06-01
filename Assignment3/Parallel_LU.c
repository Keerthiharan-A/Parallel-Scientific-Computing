#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
double f(double x) {
    return sin(5 * x);
}

int main(int argc, char *argv[]) {
    int ngangs;
    if(argc==2){
        ngangs = atoi(argv[1]);
    }
    else{
        printf("Enter number of gangs: ");
        scanf("%d",&ngangs);
    }
    int n = 1000;
    double h = 3.0 / n;
    double a[n + 1], b[n + 1], y[n + 1], x[n + 1], c[n + 1], B[n + 1], l1[n + 1] , u1[n + 1];
    // Update a, b, c, l1 on the CPU
    b[0] = 1;
    c[0] = 2;
    u1[0] = 1;
    clock_t start_time = clock();
#pragma acc parallel loop present(f) copy(B) firstprivate(h) num_gangs(ngangs)
    for (int i = 0; i <= n; i++) {
        if (i == 0) {
            B[0] = (-5*f(0)/2.0 + 2*f(h) + f(2*h)/2.0)/h;
        } else if (i < n) {
            B[i] = 3.0*(f(h*i+h) - f(h*i-h))/h;
        } else {
            B[n] = (5*f(n*h)/2.0 - 2*f((n-1)*h) - f((n-2)*h)/2.0)/h;
        }
    }
    for (int i = 1; i < n; i++) {
        a[i] = 1;
        b[i] = 4;
        c[i] = 1;
        l1[i] = a[i] / u1[i - 1];
        u1[i] = b[i] - l1[i] * c[i - 1];
    }
    b[n] = 1;
    a[n] = 2;
    l1[n] = a[n]/(u1[n-1]);
    u1[n] = b[n] - l1[n]*c[n-1];
    y[0] = B[0];
    for(int i=1; i<=n; ++i)  y[i] = B[i] - l1[i]*y[i-1];
    x[n] = y[n]/u1[n];
    for(int i=n-1; i>=0; --i) x[i] = (y[i] - c[i]*x[i+1])/ u1[i];
    clock_t end_time = clock();
    double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    for(int i=0; i<=n; ++i) printf("%lf ", x[i]);
    printf("Execution time: %.6f seconds\n", total_time);
    return 0;
}
