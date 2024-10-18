#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIZE 7680

void MultiplyMatrixVector(double *matrix, double *vector, double *result, const int N) {
    for (int i = 0; i < N; ++i) {
        result[i] = 0;
        for (int j = 0; j < SIZE; ++j) {
            result[i] += matrix[i * SIZE + j] * vector[j];
        }
    }
}

void MultiplyVectorScalar(double *vector, double scalar, double *result, const int N) {
    for (int i = 0; i < N; ++i) {
        result[i] = vector[i] * scalar;
    }
}

double DotMultiply(const double *firstVector, const double *secondVector, const int N) {
    double result = 0;
    for (int i = 0; i < N; i++) {
        result += firstVector[i] * secondVector[i];
    }
    return result;
}

void SubVectors(const double *firstVector, const double *secondVector, double *result, const int N) {
    for (int i = 0; i < N; i++) {
        result[i] = firstVector[i] - secondVector[i];
    }
}

void SumVectors(const double *firstVector, const double *secondVector, double *result, const int N) {
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
    if (size != SIZE) {
        printf("Неправильный размер матрицы в файле.\n");
        fclose(file);
        return;
    }

    // Читаем элементы матрицы из файла
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            fscanf(file, "%lf", &matrix[i * SIZE + j]);
        }
    }

    fclose(file);
}

void SecondDataFill(double *A, double *x, double *b, const int N) {
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



double GetNorm(double *vector, const int N) {
    double result = 0;
    for (int i = 0; i < N; ++i) {
        result += vector[i] * vector[i];
    }
    return sqrt(result);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    double start, end = 0;
    int match = 0;
    int numOfProccess, currentProcess;
    MPI_Comm_rank(MPI_COMM_WORLD, &currentProcess);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProccess);
    double e = 0.00001;
    double tau = 1e-3; //потом поменять

    const int block = SIZE / numOfProccess;
    const int matrixBlock = SIZE * SIZE / numOfProccess;

    double *A;
    double *b;
    double *x = (double *) malloc(sizeof(double) * SIZE);
    double *y = (double *) malloc(sizeof(double) * SIZE);
    double *tempAx = (double *) malloc(sizeof(double) * block);
    double *tempY = (double *) malloc(sizeof(double) * block);
    double *tempAy = (double *) malloc(sizeof(double) * block);
    double *tempYtau = (double *) malloc(sizeof(double) * block);
    double *tempX = (double *) malloc(sizeof(double) * block);

    double *partB = (double *) malloc(sizeof(double) * block);
    double *partA = (double *) malloc(sizeof(double) * matrixBlock);
    double NormY, NormB;
    double numerator, denumerator;
    if (currentProcess == 0) {
        A = (double *) malloc(sizeof(double) * SIZE * SIZE);
        b = (double *) malloc(sizeof(double) * SIZE);
        SecondDataFill(A, x, b, SIZE);
        start = MPI_Wtime();
    }
    MPI_Scatter(A, matrixBlock, MPI_DOUBLE, partA, matrixBlock, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    MPI_Bcast(b, SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(b, block, MPI_DOUBLE, partB, block, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(x, SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double partNormB = DotMultiply(partB, partB, block);
    MPI_Allreduce(&partNormB, &NormB, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    NormB = sqrt(NormB);

    for (int i = 0; i < 10000; i++) {
        MultiplyMatrixVector(partA, x, tempAx, block);
        SubVectors(tempAx, partB, tempY, block);
        MPI_Allgather(tempY, block, MPI_DOUBLE, y, block, MPI_DOUBLE, MPI_COMM_WORLD);

        double partNormY = DotMultiply(tempY, tempY, block);
        MPI_Allreduce(&partNormY, &NormY, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        NormY = sqrt(NormY);

        if (NormY / NormB < e) {
            match++;
            if (match == 3) {
                if (currentProcess == 0) {
                    printf("done on %d \n", i);
                }
                break;
            }

        } else {
            match = 0;
        }

        MultiplyMatrixVector(partA, y, tempAy, block);
        double localNumerator = DotMultiply(tempY, tempAy, block);
        double localDenurametor = DotMultiply(tempAy, tempAy, block);
        double localData[2] = {localNumerator, localDenurametor};
        double globalData[2] = {0, 0};
        MPI_Reduce(                      localData, globalData, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
       MPI_Allreduce(                 localData, globalData, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//        MPI_Allreduce(&localDenurametor, &denumerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (currentProcess == 0) {
            numerator = globalData[0];
            denumerator = globalData[1];
            tau = numerator / denumerator;
        }

        MultiplyVectorScalar(tempY, tau, tempYtau, block);
        SubVectors(x + block * currentProcess, tempYtau, x + block * currentProcess, block);

    }

    if (currentProcess == 0) {
        end = MPI_Wtime();
        printf("Time taken: %lf", end - start);
        free(A);
        free(b);
    }

    free(x);
    free(y);
    free(tempAx);
    free(tempAy);
    free(tempY);
    free(tempYtau);
    free(tempX);
    MPI_Finalize();
    return 0;
}