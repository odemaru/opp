#include <iostream>
#include <mpi.h>


void PrintMatrix(double* Matrix, int row, int col) {
    int i, j; // Loop variables
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++)
            printf("%7.4f ", Matrix[i * col + j]);
        printf("\n");
    }
}

void FillA(double* pMatrix, int rowCount, int colCount) {
    int k = 0;
    for (int i = 0; i < rowCount; i++){
        for (int j = 0; j < colCount; j++){
            pMatrix[i*colCount + j] = k;
            k++;
        }
    }
}

void FillB(double* pMatrix, int rowCount, int colCount) {
    int k = 0;
    for (int i = 0; i < rowCount; i++){
        for (int j = 0; j < colCount; j++){
            pMatrix[i*colCount + j] = k;
            k++;
        }
    }
}


void SetZero(double* pMatrix, int rowCount, int colCount) {
    for (int i = 0; i < rowCount; i++){
        for (int j = 0; j < colCount; j++){
            pMatrix[i*colCount + j] = 0;
        }
    }
}


void MultiplyMatrices(double* localA, double* localB, double* localC, int n1, int n2, int n3) {
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++)
            for (int k = 0; k < n2; k++)
                localC[i * n3 + j] += localA[i * n2 + k] * localB[k * n3 + j];
    }
}




int n1 = 4;
int n2 = 4;
int n3 = 4;



void InitializeProcess(double* &localA, double* &localB, double* &localC, double* &localBlockA, double* &localBlockB, double* &localBlockC, int &sizeBlockA, int &sizeBlockB, int procRank, int p1, int p2) {

    sizeBlockA = n1 / p1;
    sizeBlockB = n3 / p2;

    localBlockA = new double [n2 * sizeBlockA];
    localBlockB = new double [n2 * sizeBlockB];
    localBlockC = new double [sizeBlockA * sizeBlockB];
    SetZero(localBlockC, sizeBlockA, sizeBlockB);
    if (procRank == 0){
        localA = new double [n1 * n2];
        localB = new double [n2 * n3];
        localC = new double [n1 * n3];
        FillA(localA, n1, n2);
        FillB(localB, n2, n3);
        SetZero(localC, n1, n3);
    }
}




int main(int argc, char* argv[]) {
    if(argc < 3) {
        printf("./lab3 p1 p2\n");
        return 1;
    }

    int p1 = atoi(argv[1]);
    int p2 = atoi(argv[2]);


    double* A = NULL;
    double* B = NULL;
    double* C = NULL;

    int gridCoords[2]; //будущие координаты решетки
    int procNum = 0;
    int procRank = 0;
    MPI_Comm Grid;    // коммуникатор решетки
    MPI_Comm row;
    MPI_Comm col;

    double *blockA = NULL;
    int sizeBlockA = 0;
    double *blockB = NULL;
    int sizeBlockB = 0;
    double *blockC = NULL;

    double start, end = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    if ((n1 % p1 != 0) || (n3 % p2 != 0)) {
        if (procRank == 0) {
            printf("invalid grid size\n");
            return 1;
        }
    }


    //Cоздадим решетку
    int dimSize[2]; //количество процессов в каждом "измерении"
    dimSize[0] = p1;
    dimSize[1] = p2;
    int period[2]; //Задам переодичность
    period[0] = 0;
    period[1] = 0;
    int reorder = 0; // нельзя перенумеровывать
    int dimNum = 2; // колво измерений - размерность решетки
    MPI_Cart_create(MPI_COMM_WORLD, dimNum, dimSize, period, 0, &Grid); // создали декартову решетку
    MPI_Cart_coords(Grid, procRank, dimNum, gridCoords); //определение координат по номеру процесса

    int subDim[2]; //для создания подрешеток
    subDim[0] = 0;
    subDim[1] = 1; //если 1 значит будет в подрешетке
    MPI_Cart_sub(Grid, subDim, &col);
    subDim[0] = 1;
    subDim[1] = 0;
    MPI_Cart_sub(Grid, subDim, &row);
//    int rowRank, rowSize, colRank, colSize;
//    MPI_Comm_rank(col, &rowRank);
//    MPI_Comm_size(col, &rowSize);
//    MPI_Comm_rank(row, &colRank);
//    MPI_Comm_size(row, &colSize);

    //Выделяем память
    InitializeProcess(A, B, C, blockA, blockB, blockC, sizeBlockA, sizeBlockB, procRank, p1, p2);


    if (procRank == 0) {
        start = MPI_Wtime();
        printf("start time...\n");
    }
    //Раскидываем
    if (gridCoords[1] == 0) {
        MPI_Scatter(A, sizeBlockA * n2, MPI_DOUBLE, blockA,sizeBlockA * n2, MPI_DOUBLE, 0, row);
//        printf("scatter at %d \n", procRank);
//        PrintMatrix(blockA, n2, sizeBlockA);

    }

    MPI_Bcast(blockA, sizeBlockA * n2, MPI_DOUBLE, 0, col);
//    printf("bcast at %d \n", procRank);
//    PrintMatrix(blockA, n2, sizeBlockA);
//    printf("proc %d rowRank %d colRank %d", procRank, rowRank, colRank);
    MPI_Datatype colVector, colType;

    MPI_Type_vector(n2, sizeBlockB, n3, MPI_DOUBLE, &colVector); // создаем новый тип из блоков указанного размера (колво блоков, колво элементов в блоке, шаг, старый тип, куда)
    MPI_Type_commit(&colVector);
    MPI_Type_create_resized(colVector, 0, sizeBlockB * sizeof(double), &colType); // (старый, откуда, экстент, куда)
    MPI_Type_commit(&colType);

    if (gridCoords[0] == 0) {
        MPI_Scatter(B, 1, colType, blockB, n2 * sizeBlockB, MPI_DOUBLE, 0, col); //раскидываем B для строчек
//        printf("scatter at %d\n", procRank);
//        PrintMatrix(blockB, n3, sizeBlockB);
    }

    MPI_Bcast(blockB, sizeBlockB * n2, MPI_DOUBLE, 0, row);
//    printf("bcast at %d\n", procRank);
//    PrintMatrix(blockB, n3, sizeBlockB);
    //Счет
    MultiplyMatrices(blockA, blockB, blockC, sizeBlockA, n2, sizeBlockB);

    //Собираем
//    MPI_Datatype blockVector, blockType;
//    MPI_Type_vector(sizeBlockA, sizeBlockB, n3, MPI_DOUBLE, &blockVector);
//    MPI_Type_commit(&blockVector);
//
//    MPI_Type_create_resized(blockVector, 0, sizeBlockB * sizeof(double), &blockType);
//    MPI_Type_commit(&blockType);
//
//    // displs
//    int* displ =  new int[p1*p2]; //смещение относительно реквбаф
//    int* rcount =  new int[p1*p2]; //колво элементов от каждого процесса
//    int countBlock = 0;
//    int sizeBLock = sizeBlockA * sizeBlockB;
//    int count = 0;
//    int checked;
//    int j = 0;
//    printf("proc %d\n", procRank);
//    while (count < p1 * p2 * sizeBLock) {
//        checked = 0;
//        for (int i = 0; i < n3; i += sizeBlockB) {
//            displ[j] = countBlock;
//            rcount[j] = 1;
//            j++;
//            countBlock++;
//
//            checked++;
////            printf("iteration check i %d %d\n", i, displ[j-1]);
////            printf("rcount i %d\n", rcount[j-1]);
//        }
//        count += checked * sizeBLock;
//        countBlock += checked * (sizeBlockA - 1);
//
//        printf("checked after for %d\n", checked);
//        printf("while count countBlock %d %d\n", count, countBlock);
//    }
////    printf("proc rank %d\n", procRank);
////    for (int i = 0; i < sizeBLock; i++) {
////        printf("%7.4f ", blockC[i]);
////    }
////    printf("\n");
////    for (int i = 0; i < p1 * p2; i++) {
////        printf("%d ", displ[i]);
////    }
//    printf("\n");
//    MPI_Gatherv(blockC, sizeBLock, MPI_DOUBLE, C,rcount, displ, blockType, 0, MPI_COMM_WORLD);

    MPI_Datatype finalBlockC;
    MPI_Type_vector(sizeBlockA, sizeBlockB, n3, MPI_DOUBLE, &finalBlockC);
    MPI_Type_commit(&finalBlockC);

    if (procRank != 0) {
        MPI_Send(blockC, sizeBlockA * sizeBlockB, MPI_DOUBLE, 0, 1337, MPI_COMM_WORLD);
    } else {
        for (int i = 0; i < sizeBlockA; ++i) {
            for (int j = 0; j < sizeBlockB; ++j) {
                C[i * n3 + j] = blockC[sizeBlockB * i + j];
            }
        }
        for (int i = 0; i < p1; ++i) {
            for (int j = 0; j < p2; ++j) {
                if (i == 0 && j == 0) {
                    continue;
                }
                MPI_Recv(C + (i * sizeBlockA * n3 + j * sizeBlockB), 1, finalBlockC, i * p2 + j, 1337, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }



    if (procRank == 0){

        end = MPI_Wtime();
//        printf("matrix C \n");
//        for (int i = 0; i < n1* n3; i++) {
//            printf("%7.4f ", C[i]);
//        }
//        printf("\n");
//        PrintMatrix(C, n1, n3);
        printf("%lf seconds\n",end - start);
    }

    if (procRank == 0){
        delete [] A;
        delete [] B;
        delete [] C;
    }
    delete [] blockA;
    delete [] blockB;
    delete [] blockC;
    MPI_Finalize();
    return 0;
}