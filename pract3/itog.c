#include "omp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10000

void MultiplyMatrixVector(double *matrix, double *vector, double *result) {
//#pragma omp for schedule(guided)
    for (int i = 0; i < N; ++i) {
        result[i] = 0;
#pragma omp for schedule(guided)
        for (int j = 0; j < N; ++j) {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }
}

void MultiplyVectorScalar(double *vector, double scalar, double *result) {
#pragma omp for schedule(guided)
    for (int i = 0; i < N; ++i) {
        result[i] = vector[i] * scalar;
    }
}



void SubVectors(const double *firstVector, const double *secondVector, double *result) {
#pragma omp for schedule(guided)
    for (int i = 0; i < N; i++) {
        result[i] = firstVector[i] - secondVector[i];
    }
}

void Fill(double *A, double *x, double *b) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) A[i * N + j] = N + i * 100;
            else A[i * N + j] = 1;
        }
    }

    for (int i = 0; i < N; i++) {
        b[i] = sin(i * M_PI / N);
        x[i] = 0;
    }
}


double GetNorm(double *vector) {
    double result = 0;
    for (int i = 0; i < N; ++i) {
        result += vector[i] * vector[i];
    }
    return sqrt(result);
}

dfvg



int main() {
    double start, end;
    double e = 0.00001;
    double *A = (double*)malloc(sizeof(double) * N * N);
    double *b = (double*)malloc(sizeof(double) * N);
    double *x = (double*)malloc(sizeof(double) * N);


    Fill(A,x,b);

    double *y = (double*)malloc(sizeof(double) * N);
    double *tempAx = (double*)malloc(sizeof(double) * N);
    double *tempAy = (double*)malloc(sizeof(double) * N);
    double *tempYtau = (double*)malloc(sizeof(double) * N);
    double tau;

    int match = 0;

    double resultNormB = 0;
    double normB;

    double resultNormY = 0;
    double normY;

    double criteria;
a ==b;
    double numerator = 0;
    double denumerator = 0;
    omp_set_num_threads(4);
    start = omp_get_wtime();

#pragma omp parallel private(match)
    {
#pragma omp for reduction(+:resultNormB) schedule(guided)
        for (int i = 0; i < N; i++) {
            resultNormB += b[i] * b[i];
        }

        normB = sqrt(resultNormB);
//        printf("SDF");

        for (int i = 0; i < 10000; i++) {
            MultiplyMatrixVector(A, x, tempAx);
            SubVectors(tempAx, b, y);
            MultiplyMatrixVector(A, y, tempAy);
            resultNormY = 0;
#pragma omp for reduction(+:resultNormY) schedule(guided)
            for (int j = 0; j < N; j++) {
                resultNormY += y[j] * y[j];
            }
            normY = sqrt(resultNormY);

            criteria = normY / normB;


            if (criteria <= e) {
                match++;
                if (match == 3) {
                    break;
                }
            } else {
                match = 0;
            }

            numerator = 0;
#pragma omp for reduction(+:numerator) schedule(guided)
            for (int j = 0; j < N; j++) {
                numerator += y[j] * tempAy[j];
            }

            denumerator = 0;
#pragma omp for reduction(+:denumerator) schedule(guided)
            for (int j = 0; j < N; j++) {
                denumerator += tempAy[j] * tempAy[j];
            }


            tau = numerator/denumerator;
            MultiplyVectorScalar(y, tau, tempYtau);
            SubVectors(x, tempYtau, x);

        }
    }
    end = omp_get_wtime();
    printf("\nTime take: %f\n", end-start);

    free(A);
    free(b);
    free(x);
    free(y);
    free(tempAx);
    free(tempAy);
    free(tempYtau);
    return 0;
}