#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10000

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

void MultiplyMatrixFullVector(double *matrix, double *vector, double *result, int size, int sizeBlock) {
    for (int i = 0; i < sizeBlock; ++i) {
        result[i] = 0;
        for (int j = 0; j < size; ++j) {
            result[i] += matrix[i * size + j] * vector[j];
        }
    }
}

double DotMultiply(double *firstVector, double *secondVector, int size) {
    double result = 0;
    for (int i = 0; i < size; ++i) {
        result += firstVector[i] * secondVector[i];
    }
    return result;
}

void SubVectors(const double *firstVector, const double *secondVector, double *result, int size) {
    for (int i = 0; i < size; i++) {
        result[i] = firstVector[i] - secondVector[i];
    }
}

void MultiplyVectorScalar(double *vector, double scalar, double *result, int size) {
    for (int i = 0; i < size; ++i) {
        result[i] = vector[i] * scalar;
    }

}

int main(int argc, char**argv) {
    MPI_Init(&argc, &argv);
    double start, end;
    int rank, numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double e = 0.00001;
    double tau;
    double *x = (double*)malloc(sizeof(double) * N);
    double *y = (double*)malloc(sizeof(double) * N);
    double *mainA;
    double *mainB = (double*)malloc(sizeof(double) * N);

    if (rank == 0) {
        mainA = (double*)malloc(sizeof(double) * N * N);
        Fill(mainA, x, mainB);
        start = MPI_Wtime();
    }

    //Раздача
    int matrixBlock = N * N / numProc;
    int block = N / numProc;

    double *blockA = (double*)malloc(sizeof(double) * matrixBlock); //разрезаем и раздаем всем процессам А
    MPI_Scatter(mainA, matrixBlock, MPI_DOUBLE, blockA, matrixBlock, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //разрезаем и раздаем Б
//    double* blockB = (double*)malloc(sizeof(double) * block);
//    MPI_Scatter(mainB, block, MPI_DOUBLE, blockB, block, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(mainB, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //все процессы будут иметь ПОЛНЫЙ вектор Х
    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Темпы
    double *blockAx = (double*)malloc(sizeof(double) * block);
    double *blockY = (double*)malloc(sizeof(double) * block);
    double *blockAy = (double*)malloc(sizeof(double) * block);
    double *blockYtau = (double*)malloc(sizeof(double) * block);

    //Критерий normY/normB < epsilon
    double normY = 0;
    double normB = 0;
//    double blockNormB = DotMultiply(blockB, blockB, block);
//    MPI_Allreduce(&blockNormB, &normB, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    normB = DotMultiply(mainB, mainB, N);
    normB = sqrt(normB);
    int match = 0;
    double numerator, denumerator;
    for (int i = 0; i < 10000; i++) {
        MultiplyMatrixFullVector(blockA, x, blockAx, N, block); //считаем Ах
        SubVectors(blockAx, &mainB[rank * block], blockY, block); //cчитаем Ax-b

        //Критерий в будущем
        double blockNormY = DotMultiply(blockY, blockY, block);
        MPI_Allreduce(&blockNormY, &normY, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        normY = sqrt(normY);// со всех процессов складываем полную норму Y

        //Критерий
        if (normY / normB < e) {
            match++;
            if (match == 3) {
                if (rank == 0) {
                    printf("Done at %d\n", i);
                }
                break;
            }
        } else {
            match = 0;
        }
        // Собираем вектор Y со всех процессов (блок 0 процесса в 0 блок буффера)
        MPI_Allgather(blockY, block, MPI_DOUBLE, y, block, MPI_DOUBLE, MPI_COMM_WORLD);
        // Считаем Ay
        MultiplyMatrixFullVector(blockA, y, blockAy, N, block);
        double blockNumerator = DotMultiply(blockY, blockAy, block); //числитель Тау
        double blockDenumerator = DotMultiply(blockAy, blockAy, block); //знаменатель Тау
        double localData[2] = {blockNumerator, blockDenumerator};
        double fullData[2] = {0, 0};
        MPI_reduce(localData, fullData, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        numerator = fullData[0];
        denumerator = fullData[1];
        tau = numerator / denumerator;

        MultiplyVectorScalar(blockY, tau, blockYtau, block);
        SubVectors(x + block * rank, blockYtau, x + block * rank, block);

    }

    if (rank == 0) {
        if (match < 3) {
            printf("ERROR");
        } else {
            end = MPI_Wtime();
            printf("Time taken: %lf\n", end - start);
        }

        free(mainA);
        free(mainB);
    }

    free(x);
    free(y);

    free(blockY);
    free(blockAy);
    free(blockA);
//    free(blockB);
    free(blockYtau);
    free(blockAx);
    MPI_Finalize();
    return 0;
}