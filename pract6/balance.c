#include <pthread.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SEND_TASKS 2
#define SEND_TASKS_COUNT 3

#define L 70
#define GLOBAL_NUM_ITERATIONS 100
#define NUM_OF_TASKS 50
#define MIN_TASKS_COUNT 2

#define DONE -1
#define NO_TASKS -2

int procNum;
int procRank;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

int *tasks;

double disbalance = 0;
int finishStatus = 0;

int numOfIterationTasks, extraTaskCount, numOfDoneTasks;
double result = 0;

void MainCalculation(const int *tasks) {
    for(int i = 0; i < numOfIterationTasks; i++) {
        pthread_mutex_lock(&mutex);
        int repeatNum = tasks[i];
        pthread_mutex_unlock(&mutex);

        numOfDoneTasks++;
        for(int j = 0; j < repeatNum; j++) {
            result += sin(j);
        }
    }
    numOfIterationTasks = 0;
}

void Fill(int *tasks, int numTasks, int iterCounter) {
    for (int i = 0; i < numTasks; i++)
        tasks[i] = abs(50 - i % 100) * abs(procRank - (iterCounter % procNum)) * L;
}

void* CalculationByThread(void *me) {
    double duration, minDuration, maxDuration;

    tasks = (int*)malloc(sizeof(double) * NUM_OF_TASKS);
    printf("START CALC");
    double startThread, endThread;
    for(int i = 0; i < GLOBAL_NUM_ITERATIONS; i++) {
        startThread = MPI_Wtime();
        printf("ITER %d", i);
        MPI_Barrier(MPI_COMM_WORLD);

        Fill(tasks, NUM_OF_TASKS, i);

        numOfIterationTasks = NUM_OF_TASKS;

        numOfDoneTasks = 0;
        extraTaskCount = 0;

        MainCalculation(tasks);

        int response;
        for (int process = 0; process < procNum; ++process) {
            if (process != procRank){
                MPI_Send(&procRank, 1, MPI_INT, process, 0, MPI_COMM_WORLD);
                MPI_Recv(&response, 1, MPI_INT, process, SEND_TASKS_COUNT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (response != NO_TASKS){
                    extraTaskCount = response;
                    memset(tasks, 0, NUM_OF_TASKS);

                    MPI_Recv(tasks, extraTaskCount, MPI_INT, process, SEND_TASKS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    pthread_mutex_lock(&mutex);
                    numOfIterationTasks = extraTaskCount;
                    pthread_mutex_unlock(&mutex);

                    MainCalculation(tasks);
                }
            }
        }
        endThread = MPI_Wtime();
        duration = endThread - startThread;

        MPI_Allreduce(&duration, &maxDuration, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&duration, &minDuration, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        disbalance += (maxDuration - minDuration) / maxDuration;

    }

    pthread_mutex_lock(&mutex);
    finishStatus = 1;
    pthread_mutex_unlock(&mutex);

    int signature = DONE;
    MPI_Send(&signature, 1, MPI_INT, procRank, 0, MPI_COMM_WORLD);

    free(tasks);
    pthread_exit(NULL);
}

void* receiverThread(void *me) {
    int requestRank, sendTaskNum, message;

    MPI_Status status;
    MPI_Barrier(MPI_COMM_WORLD);

    while (finishStatus != 1){
        MPI_Recv(&message, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

        if (message == DONE)
            printf("at proc %d done all tasks\n", procRank);

        requestRank = message;

        pthread_mutex_lock(&mutex);
        if (numOfIterationTasks >= MIN_TASKS_COUNT){
            sendTaskNum = numOfIterationTasks / (procNum * 2);
            numOfIterationTasks = numOfIterationTasks / (procNum * 2);
            MPI_Send(&sendTaskNum, 1, MPI_INT, requestRank, SEND_TASKS_COUNT, MPI_COMM_WORLD);
            MPI_Send(&tasks[NUM_OF_TASKS - sendTaskNum], sendTaskNum, MPI_INT, requestRank, SEND_TASKS, MPI_COMM_WORLD);
        } else {
            sendTaskNum = NO_TASKS;
            MPI_Send(&sendTaskNum, 1, MPI_INT, requestRank, SEND_TASKS_COUNT, MPI_COMM_WORLD);
        }
        pthread_mutex_unlock(&mutex);

    }
    pthread_exit(NULL);
}


int main(int argc, char* argv[]) {
    int threadSupport;

    pthread_t thr[2];
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &threadSupport);

    if (threadSupport != MPI_THREAD_MULTIPLE) {
        MPI_Finalize();
        return -1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);

    pthread_mutex_init(&mutex, NULL);
    pthread_attr_t ThreadAttributes;

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    printf("INIT");

    pthread_attr_init(&ThreadAttributes);
    pthread_attr_setdetachstate(&ThreadAttributes, PTHREAD_CREATE_JOINABLE);

    pthread_create(&thr[0], &ThreadAttributes, receiverThread, NULL);
    pthread_create(&thr[1], &ThreadAttributes, CalculationByThread, NULL);
    printf("JOIN");

    pthread_join(thr[0], NULL);
    pthread_join(thr[1], NULL);

    pthread_attr_destroy(&ThreadAttributes);
    pthread_mutex_destroy(&mutex);

    double end = MPI_Wtime();
    if (procRank == 0){
        printf("Time taken: %f\n" , end - start);
        printf("Total disbalance: %f\n", disbalance / GLOBAL_NUM_ITERATIONS * 100);
    }

    MPI_Finalize();
    return 0;
}
