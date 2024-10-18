#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10000

void MultiplyMatrixVector(double *matrix, double *vector, double *result) {
    for (int i = 0; i < N; ++i) {
        result[i] = 0;
        for (int j = 0; j < N; ++j) {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }
}

void MultiplyVectorScalar(double *vector, double scalar, double *result) {
    for (int i = 0; i < N; ++i) {
        result[i] = vector[i] * scalar;
    }
}

double DotMultiply(const double *firstVector, const double *secondVector) {
    double result = 0;
    for (int i = 0; i < N; i++) {
        result += firstVector[i] * secondVector[i];
    }
    return result;
}

void SubVectors(const double *firstVector, const double *secondVector, double *result) {
    for (int i = 0; i < N; i++) {
        result[i] = firstVector[i] - secondVector[i];
    }
}

void SumVectors(const double *firstVector, const double *secondVector, double *result) {
    for (int i = 0; i < N; i++) {
        result[i] = firstVector[i] + secondVector[i];
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





int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    const double start = MPI_Wtime();
    int numOfProccess, currentProcess;
    MPI_Comm_rank(MPI_COMM_WORLD, &currentProcess);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProccess);
    double e = 0.00001;
    double *A = (double*)malloc(sizeof(double) * N * N);
    double *b = (double*)malloc(sizeof(double) * N);
    double *x = (double*)malloc(sizeof(double) * N);


    SecondDataFill(A,x,b);

    double *y = (double*)malloc(sizeof(double) * N);
    double *tempAx = (double*)malloc(sizeof(double) * N);
    double *tempAy = (double*)malloc(sizeof(double) * N);
    double *tempYtau = (double*)malloc(sizeof(double) * N);

    int match = 0;

    for (int i = 0; i < 10000; i++) {
        MultiplyMatrixVector(A, x, tempAx);
        SubVectors(tempAx, b, y);
        MultiplyMatrixVector(A, y, tempAy);

        double criteria = GetNorm(y) / GetNorm(b);


        if (criteria < e) {
            match++;
            if (match == 3) {
                printf("done\n");
                break;
            }
        } else {
            match = 0;
        }

        double numerator = DotMultiply(y, tempAy);
        double denominator = DotMultiply(tempAy, tempAy);
        double tau = numerator/denominator;


        MultiplyVectorScalar(y, tau, tempYtau);
        SubVectors(x, tempYtau, x);

    }

    if (match < 3) {
        printf("wtf");
    } else {
        const double end = MPI_Wtime();
        printf("Time taken: %lf", end - start);
    }

    free(A);
    free(b);
    free(x);
    free(y);
    free(tempAx);
    free(tempAy);
    free(tempYtau);
    MPI_Finalize();
    return 0;
}