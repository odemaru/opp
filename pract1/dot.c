#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define FIRST_TAG 0
#define SECOND_TAG 1
#define FINAL_TAG 2


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
    MPI_Status status;

    //0 процесс
    if (curProcess == 0) {
        firstVector = CreateVector(length);
        secondVector = CreateVector(length);

        unsigned long sum, partSum = 0;
        unsigned long part = length / numOfProcess;
        unsigned long shift = length % numOfProcess;
        start = MPI_Wtime();
        for (int i = 1; i < numOfProcess; ++i){
            MPI_Send(firstVector + part * i + shift, part, MPI_INT, i, FIRST_TAG, MPI_COMM_WORLD);
            // сcылка на 1 элем (где данные лежат), колво элементов, тип данных, кому отправляем номер процесса, тэг сообщение, коммуникатор
            MPI_Send(secondVector, length, MPI_INT, i, SECOND_TAG, MPI_COMM_WORLD);
        }

        SumOfVectors(firstVector, secondVector, shift + part, length, &sum);

        for (int i = 1; i < numOfProcess; ++i) {
            MPI_Recv(&partSum, 1, MPI_UNSIGNED_LONG, i, 2, MPI_COMM_WORLD, &status);
            //ссылка куда сохранять, макс колво, тип данных, номер процесса который отправил, тэг сообщения, коммутатор, статус
            sum += partSum;
        }
        end = MPI_Wtime();

        printf("Sum: %ld\nTime: %lf seconds\n", sum, end - start);
    } else { //остальные процессы
        unsigned long sum = 0;
        int part = length / numOfProcess;

        firstVector = (int*)malloc(sizeof(int) * part);
        secondVector = (int*)malloc(sizeof(int) * length);

        MPI_Recv(firstVector, part, MPI_INT, 0, FIRST_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(secondVector, length, MPI_INT, 0, SECOND_TAG, MPI_COMM_WORLD, &status);

        SumOfVectors(firstVector, secondVector, part, length, &sum);
        MPI_Send(&sum, 1, MPI_UNSIGNED_LONG, 0, FINAL_TAG, MPI_COMM_WORLD);
    }

    free(firstVector);
    free(secondVector);

    MPI_Finalize();
    return 0;
}