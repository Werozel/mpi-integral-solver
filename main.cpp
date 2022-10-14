#include <iostream>
#include "mpi.h"

int main(int argc, char *argv[]) {
    int size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (!rank) {
        printf("Hello world\n");
        fflush(stdin);
    } else {
        printf("Hello from %d\n", rank);
        fflush(stdin);
    }

    MPI_Finalize();

    return 0;
}
