#include <iostream>
#include <mpi.h>
#include <cmath>

double a1 = 0, a2 = 0, a3 = 0, b1 = 1, b2 = 1, b3 = 1;
double area = (b1 - a1) * (b2 - a2) * (b3 - a3);

const int NUMBER_OF_POINTS = 1000000;
const int arr_size = NUMBER_OF_POINTS * 3;
int n_proc, rank;


double f(double x, double y) {
    // TODO: check if inside area
    return sqrt(x * x + y * y);
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void master_routine() {
    srand(time(nullptr));

    double *points = (double *) malloc(arr_size * sizeof(double) + MPI_BSEND_OVERHEAD + 7);
    for (int i = 0; i < arr_size; i += 3) {
        double x = fRand(a1, b1);
        double y = fRand(a2, b2);
        double z = fRand(a3, b3);

        points[i] = x;
        points[i + 1] = y;
        points[i + 2] = z;
    }

    int bunch_size = NUMBER_OF_POINTS / (n_proc - 1) * 3;
    for (int i = 1; i < n_proc; i++) {
        int offset = (i - 1) * bunch_size;
        MPI_Send(&points[offset], bunch_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }

    double sum = 0;
    for (int i = 1; i < n_proc; i++) {
        double local_sum;
        MPI_Recv(&local_sum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        sum += local_sum;
    }

    double result = area * sum / NUMBER_OF_POINTS;
    std::cout << result << std::endl;
}

void worker_routine() {
    int bunch_size = NUMBER_OF_POINTS / (n_proc - 1) * 3;
    double *points = (double *) malloc(bunch_size * sizeof(double) + MPI_BSEND_OVERHEAD + 7);

    MPI_Recv(points, bunch_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    double local_sum = 0;
    for (int i = 0; i < bunch_size; i += 3) {
        local_sum += f(points[i], points[i + 1]);
    }

    MPI_Send(&local_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}

int main(int argc, char *argv[]) {
    int err = MPI_Init(&argc, &argv);

    if (err != MPI_SUCCESS)
    {
        fprintf(stderr, "Error while starting! \n");
        MPI_Abort(MPI_COMM_WORLD, err);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

    if (!rank) {
        master_routine();
    } else {
        worker_routine();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
