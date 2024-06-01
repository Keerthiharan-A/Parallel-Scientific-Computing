#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TYPE        float
#define N           100
#define SMALLVALUE  0.001

void printMat(TYPE **matrix, int size) {

    #pragma acc serial loop present(matrix[0:size][0:size])
    for (int i = 0; i < size; i++) {
        for (int j = 0; j <= i; j++)
            printf("%8.6f ", matrix[i][j]);
        printf("\n");
    }
    return;
}

void init(TYPE **mat, int size) {

#pragma acc parallel loop present(mat[0:size][0:size]) collapse(2)
    for (int ii = 0; ii < size; ++ii)
        for (int jj = 0; jj < size; ++jj)
            mat[ii][jj] = (ii + jj) / (float)size / size;

#pragma acc parallel loop present(mat[0:size][0:size])
    for (int ii = 0; ii < size; ++ii)
        mat[ii][ii] = 1.0;
    return;
}

void cholesky(TYPE **a, int size) {
    for (int ii = 0; ii < size; ++ii) {
        float sum = a[ii][ii];
#pragma acc parallel loop present(a[0:size][0:size]) firstprivate(ii) reduction(+:sum)
        for (int kk = 0; kk < ii; ++kk)
            sum += -a[ii][kk] * a[ii][kk];
        a[ii][ii] = sqrt(sum);
#pragma acc parallel loop present(a[0:size][0:size], SMALLVALUE) firstprivate(ii)
        for (int jj = ii + 1; jj < size; ++jj) {
            for (int kk = 0; kk < ii; ++kk)
                a[jj][ii] += -a[jj][kk] * a[ii][kk];
            a[jj][ii] /= (a[ii][ii] > SMALLVALUE ? a[ii][ii] : 1);
        }
    }
    return;
}

int main() {

    TYPE **a = (TYPE **)malloc(N * sizeof(TYPE *));
    for (int i = 0; i < N; i++)
        a[i] = (TYPE *)malloc(N * sizeof(TYPE));

    clock_t start_time = clock();
#pragma acc data copyin(a[0:N][0:N], SMALLVALUE)
    {
        init(a, N);
        cholesky(a, N);
        printMat(a, N);
    }
    clock_t end_time = clock();
    double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %.6f seconds\n", total_time);

    // Free dynamically allocated memory
    for (int i = 0; i < N; i++)
        free(a[i]);
    free(a);

    return 0;
}
