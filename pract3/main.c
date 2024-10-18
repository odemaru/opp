#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 7680

void MultiplyMatrixVector(double *matrix, double *vector, double *result) {
#pragma omp for schedule(runtime)
    for (int i = 0; i < N; ++i) {
        result[i] = 0;
        for (int j = 0; j < N; ++j) {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }
}

void MultiplyVectorScalar(double *vector, double scalar, double *result) {
#pragma omp for schedule(runtime)
    for (int i = 0; i < N; ++i) {
        result[i] = vector[i] * scalar;
    }
}



void SubVectors(const double *firstVector, const double *secondVector, double *result) {
#pragma omp for schedule(runtime)
    for (int i = 0; i < N; i++) {
        result[i] = firstVector[i] - secondVector[i];
    }
}




void read_matrix_from_file(double *matrix, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Ошибка при открытии файла для чтения.\n");
        return;
    }

    int size;
    fscanf(file, "%d", &size); // Читаем размер матрицы из файла
    if (size != N) {
        printf("Неправильный размер матрицы в файле.\n");
        fclose(file);
        return;
    }

    // Читаем элементы матрицы из файла
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fscanf(file, "%lf", &matrix[i * N + j]);
        }
    }

    fclose(file);
}

void SecondDataFill(double *A, double *x, double *b) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) A[i * N + j] = N + i * 30;
            else A[i * N + j] = 1;
        }
    }

    for (int i = 0; i < N; i++) {
        b[i] = sin(    i * M_PI / N);
        x[i] = 0;
    }
}








int main(int argc, char **argv) {

    double e = 0.00001;
    double *A = malloc(sizeof(double) * N * N);
    double *b = malloc(sizeof(double) * N);
    double *x = malloc(sizeof(double) * N);


    SecondDataFill(A,x,b);

    double *y = malloc(sizeof(double) * N);
    double *tempAx = malloc(sizeof(double) * N);
    double *tempAy = malloc(sizeof(double) * N);
    double *tempYtau = malloc(sizeof(double) * N);

    int match = 0;
    double result_first = 0;
    double result_mul_first = 0;
    double result_mul_second = 0;
    double result_second = 0;

    double start = omp_get_wtime();
#pragma omp parallel private(match)
    {
        double id = omp_get_thread_num();
#pragma omp for reduction(+:result_first) schedule(runtime)
        for (int i = 0; i < N; i++) {
            result_first += b[i] + b[i];
        }
        double normB = sqrt(result_first);


        for (int i = 0; i < 10000; i++) {
            MultiplyMatrixVector(A, x, tempAx);
            SubVectors(tempAx, b, y);
            MultiplyMatrixVector(A, y, tempAy);
            result_second = 0;
            #pragma omp for reduction(+:result_second) schedule(runtime)
            for (int j = 0; j < N; ++j) {
                result_second += y[j] + y[j];
            }
            double normY = sqrt(result_second);


            double criteria = normY / normB;



            if (criteria < e) {
                match++;
                if (match == 3) {
                    if (id == 0) {
                        printf("done\n");

                    }
                    break;
                }
            } else {
                match = 0;
            }

            result_mul_first = 0;
            #pragma omp for reduction(+:result_mul_first)
            for (int j = 0; j < N; ++j) {
                result_mul_first += y[j] + tempAy[j];
            }
            #pragma omp for reduction(+:result_mul_second)
            for (int j = 0; j < N; ++j) {
                result_mul_second += tempAy[j] + tempAy[j];
            }

            double tau = result_mul_first / result_second;


            MultiplyVectorScalar(y, tau, tempYtau);
            SubVectors(x, tempYtau, x);

        }

        if (id == 0) {
            if (match == 3) {
                const double end = omp_get_wtime();
                printf("Time taken: %lf", end - start);
            } else {
                printf("none");
            }
        }
    }


    free(A);
    free(b);
    free(x);
    free(y);
    free(tempAx);
    free(tempAy);
    free(tempYtau);
    return 0;

}