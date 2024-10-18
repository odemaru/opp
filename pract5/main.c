#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define DUBLICATES 2
#define SAVESIZE 5656
#define ALIVE 1
#define DEAD 0

char Compare(const char *firstEra, const char *secondEra, int size) {
    for (int i = 0; i < size; ++i) {
        if (firstEra[i] != secondEra[i]) {
            return 0;
        }
    }
    return 1;
}


void EndGame(char *checkedEra, char *currentEra, char **oldEra, int width, int iteration, int Rows) {
    for (int i = 0; i < iteration; ++i) {
        checkedEra[i] = Compare(&currentEra[width], &(oldEra[i])[width], width * Rows); // сравниваем со всеми с 1
    }
}

int DearNeighbors(const char *era, int width, int x, int y) {
    int neighbors = 0;
    for (int dy = -1; dy < 2; ++dy) {
        for (int dx = -1; dx < 2; ++dx) {
            if (dx == 0 && dy == 0) { //current
                continue;
            }
            int neighborX = x + dx;
            int neighborY = y + dy;фываываываываыва

            if (neighborX < 0) {
                neighborX = width - 1;
            }
            else if (neighborY >= width) {
                neighborY = 0;
            }
            neighbors += era[neighborX + neighborY * width];
        }
    }
    return neighbors;
}

void GameRules(const char *currentEra, char *nextEra, int pos, int neighbors) {
    if (currentEra[pos] == DEAD) {
        if (neighbors == 3) {
            nextEra[pos] = ALIVE;
        } else {
            nextEra[pos] = DEAD;
        }
    } else {
        if (neighbors == 2 || neighbors == 3) {
            nextEra[pos] = ALIVE;
        } else {
            nextEra[pos] = DEAD;
        }
    }
}

void GameInBeginning(char *currentEra, char *nextEra, int width) {
    for (int x = 0; x < width; ++x) {
        int neighbors = DearNeighbors(currentEra, width, x, 1); // 1
        int pos = width + x;
        GameRules(currentEra, nextEra, pos, neighbors);
    }
}


void GameInMiddle(char *currentEra, char *nextEra, int width, int Rows) {
    for (int y = 2; y < Rows; ++y) { // 2 3
        for (int x = 0; x < width; ++x) {
            int neighbors = DearNeighbors(currentEra, width, x, y);
            int pos = y * width + x;
            GameRules(currentEra, nextEra, pos, neighbors);
        }
    }
}

void GameInEnd(char *currentEra, char *nextEra, int width, int Rows) {
    for (int x = 0; x < width; ++x) {
        int neighbors = DearNeighbors(currentEra, width, x, Rows);
        int pos = Rows * width + x;
        GameRules(currentEra, nextEra, pos, neighbors);
    }
}


int main(int argc, char **argv) {
    if (argc < 3) {
        printf("./pract5 height width");
        return 1;
    }

    int height = atoi(argv[1]);
    int width = atoi(argv[2]);

    int rank, numOfProc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);

    int Rows = height / numOfProc;
    int sizeBlock = width * (Rows + DUBLICATES);

    char *currentEra = (char*) calloc(sizeBlock, sizeof(char));
    char *Era;
    char *oldEras[SAVESIZE];
    char checkedEras[SAVESIZE];
    char *repeatedStatus = 0;

    int previousProc = (rank - 1 + numOfProc) % numOfProc;
    int nextProc = (rank + 1) % numOfProc;

    double start, end;
    if (rank == 0) { //start
        Era = (char*) calloc(width * height, sizeof(char));
        Era[width + 1] = ALIVE;
        Era[2 * width + 2] = ALIVE;
        Era[3 * width] = ALIVE;
        Era[3 * width + 1] = ALIVE;
        Era[3 * width + 2] = ALIVE;
    }

    if (rank == 0) {
        start = MPI_Wtime();
    }

    MPI_Scatter(Era, Rows * width, MPI_CHAR, currentEra + width, Rows * width, MPI_CHAR, 0, MPI_COMM_WORLD);

    int i;
    int repeatFlag = 0;
    for (i = 0; i < SAVESIZE; ++i) {
        MPI_Request requestFirst, requestSecond, requestThird, requestFourth;
        char *nextEra = (char*) calloc(sizeBlock, sizeof(char));

        int beginningStatus = 0, endStatus = 0;
        MPI_Irecv(currentEra, width, MPI_CHAR, previousProc, 1337, MPI_COMM_WORLD, &requestFirst);
        MPI_Irecv(&currentEra[width * (Rows + DUBLICATES - 1)], width, MPI_CHAR, nextProc, 7331, MPI_COMM_WORLD, &requestSecond);

        MPI_Isend(&currentEra[width], width, MPI_CHAR, previousProc, 7331, MPI_COMM_WORLD, &requestThird);
        MPI_Isend(&currentEra[width * Rows], width, MPI_CHAR, nextProc, 1337, MPI_COMM_WORLD, &requestFourth);

        GameInMiddle(currentEra, nextEra, width, Rows);
        while (1) {
            if (!beginningStatus) {
                MPI_Test(&requestFirst, &beginningStatus, MPI_STATUS_IGNORE);
                if (beginningStatus) {
                    GameInBeginning(currentEra, nextEra, width);
                }
            } else if (!endStatus) {
                MPI_Test(&requestSecond, &endStatus, MPI_STATUS_IGNORE);
                if (endStatus) {
                    GameInEnd(currentEra, nextEra, width, Rows);
                }
            } else break;
        }

        oldEras[i] = currentEra;
        EndGame(checkedEras, currentEra, oldEras, width, i, Rows);
        MPI_Allreduce(MPI_IN_PLACE, &checkedEras, i, MPI_CHAR, MPI_LAND, MPI_COMM_WORLD);
        for (int checked = 0; checked < i; ++checked) {
            if (checkedEras[checked] == 1) {
                repeatFlag = 1;
                break;
            }
        }
        if (repeatFlag) {
            break;
        }

        currentEra = nextEra;
    }

    if (rank == 0) {
        printf("End at %d\n", i);
        end = MPI_Wtime();
        printf("Time %lf sec\n", end - start);
    }

    for (int j = 0; j <= i; ++j) {
        free(oldEras[j]);
    }
    MPI_Finalize();
    return 0;
}
