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

void SumOfVectors(int *firstVector, int *secondVector, unsigned long firstLength, unsigned long secondLength, unsigned long *sum){
    for (int i = 0; i < firstLength; ++i){
        for (int j = 0; j < secondLength; ++j){
            *sum += firstVector[i] * secondVector[j];
        }
    }
}

int main(int argc, char **argv){
    int curProcess = 0;
    int numOfProcess = 0;

    double start = 0.0;
    double end = 0.0;

    int *firstVector;
    int *secondVector;

    unsigned long length = atoi(argv[1]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &curProcess);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcess);

    unsigned long sum, partSum = 0;
    unsigned long part = length / numOfProcess;
    unsigned long shift = length % numOfProcess;
    unsigned long localSum = 0;

    int *partSecondVector = (int*)malloc(sizeof(int) * part);
    firstVector = (int*)malloc(sizeof(int)*length);

    if (curProcess == 0){
        firstVector = CreateVector(length);
        secondVector = CreateVector(length);

        start = MPI_Wtime();
        SumOfVectors(firstVector, secondVector, length, shift, &localSum);






















        \
    }

    MPI_Scatter(secondVector + shift, part, MPI_INT, partSecondVector, part, MPI_INT, 0, MPI_COMM_WORLD);// распределяет из одного члена группы по всем (нарезвает на куски)
    //ссылка что берем, колво, тип, куда кладем, колво принимаемого, тип принимаемого, процесс, коммуникатор
    MPI_Bcast(firstVector, length, MPI_INT, 0, MPI_COMM_WORLD); //транслирует из одного остальным
    //ссылка что берем, колво, тип, процесс, коммуникатор

    SumOfVectors(firstVector, partSecondVector, length, part, &localSum);

    MPI_Reduce(&localSum, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD); //собрать результаты

    if (curProcess == 0){
        end = MPI_Wtime();
        printf("Sum: %ld\nTime: %lf seconds\n", sum, end - start);

        free(firstVector);
        free(secondVector);
        free(partSecondVector);
    } else {
        free(partSecondVector);
        free(firstVector);
    }
    MPI_Finalize();
    return 0;
}