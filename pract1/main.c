#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int *CreateVector(unsigned long length){ //создаем вектор и заполняем рандомными цифрами от 1 до 9
    srand(time(NULL));
    int *vector = (int*)malloc(sizeof(int) * length);
    for (int i = 0; i < length; ++i){
        vector[i] = rand()%9 + 1;
    }

    return vector;
}

void SumOfVectors(int *firstVector, int *secondVector, unsigned long length, unsigned long *sum){
    for (int i = 0; i < length; ++i){
        for (int j = 0; j < length; ++j){
            *sum += firstVector[i] * secondVector[j];
        }
    }
}

int main(int agrc, char **argv){
    int curProcess = 0;
    int numOfProcess = 0;
    double start = 0.0;
    double end = 0.0;
    unsigned long sum = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &curProcess);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcess);

    unsigned long length = atoi(argv[1]);

    int *firstVector = CreateVector(length);
    int *secondVector = CreateVector(length);

    start = MPI_Wtime();
    SumOfVectors(firstVector, secondVector, length, &sum);
    end = MPI_Wtime();

    if (curProcess == 0){
        printf("Sum: %ld\nTime: %lf seconds", sum, end - start);
    }

    free(secondVector);
    free(firstVector);
    MPI_Finalize();
    return 0;
}