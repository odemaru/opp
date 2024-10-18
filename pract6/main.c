#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


#define L 1000
#define GLOBAL_NUM_ITERATIONS 500
#define NUM_TASKS 1500

int procNum;
int procRank;

int *tasks;

double disbalance = 0;

int numOfIterationTasks, numOfDoneTasks;
double result = 0;

void MainCalc(int *tasks) {
    for(int i = 0; i < numOfIterationTasks; i++) {
        int repeatNum = tasks[i];
        numOfDoneTasks++;
        for(int j = 0; j < repeatNum; j++) {
            result += sin(j);
        }
    }
    numOfIterationTasks = 0;
}

void FillTasks(int *tasks, int numTasks, int iterCounter) {
    for (int i = 0; i < numTasks; i++)
        tasks[i] = abs(50 - i % 100) * abs(procRank - (iterCounter % procNum)) * L;
}

void* CalculationByThread() {
    double duration, minDuration, maxDuration;
    tasks = (int*)malloc(sizeof(int) * NUM_TASKS);
    double startThread, endThread;

    for (int iteration = 0; iteration < GLOBAL_NUM_ITERATIONS; iteration++) {

        startThread = MPI_Wtime();

        MPI_Barrier(MPI_COMM_WORLD);

        FillTasks(tasks, NUM_TASKS, iteration);

        numOfIterationTasks = NUM_TASKS;
        numOfDoneTasks = 0;

        MainCalc(tasks);

        endThread = MPI_Wtime();

        duration = endThread - startThread;
        MPI_Allreduce(&duration, &maxDuration, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&duration, &minDuration, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        disbalance += (maxDuration - minDuration) / maxDuration;

    }
    free(tasks);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);

    double start = MPI_Wtime();
    CalculationByThread();
    double end = MPI_Wtime();

    if (procRank == 0){
        printf("Time taken: %f\n" , end - start);
        printf("Total disbalance: %f\n", disbalance / GLOBAL_NUM_ITERATIONS * 100);
    }

    MPI_Finalize();
    return 0;
}
р