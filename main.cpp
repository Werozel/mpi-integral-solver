#include <iostream>
#include <mpi.h>
#include <cmath>

double a1 = -1, a2 = -1, a3 = 0, b1 = 1, b2 = 1, b3 = 1;
double area = (b1 - a1) * (b2 - a2) * (b3 - a3);

const long NUMBER_OF_POINTS = 1024;
const long arr_size = NUMBER_OF_POINTS * 3;
int n_proc, rank;
double eps;
const double precise_res = M_PI / 6.0;

bool is_point_in_area(double x, double y, double z) {
    return x*x + y*y <= z*z;
}

double f(double x, double y, double z) {
    if (!is_point_in_area(x, y, z)) {
        return 0;
    }

    return sqrt(x*x + y*y);
}

double get_rand(double min, double max)
{
    double f = (double)rand() / RAND_MAX;
    return min + f * (max - min);
}

void master_routine() {
    srand(42);

    std::cout << "Precise result: " << precise_res << std::endl << std::endl;

    int exiting = 0;
    double result = 0;
    double total_sum = 0;
    double total_points_c = 0;
    double iteration = 0;
    double err;

    while (!exiting) {
        iteration++;
        double *points = (double *) malloc(arr_size * sizeof(double) + MPI_BSEND_OVERHEAD + 7);
        for (long i = 0; i < arr_size; i += 3) {
            double x = get_rand(a1, b1);
            double y = get_rand(a2, b2);
            double z = get_rand(a3, b3);

            points[i] = x;
            points[i + 1] = y;
            points[i + 2] = z;
        }

        long bunch_size = NUMBER_OF_POINTS / (n_proc - 1) * 3;
        for (long i = 1; i < n_proc; i++) {
            long offset = (i - 1) * bunch_size;
            MPI_Send(&points[offset], bunch_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }

        double local_sum = 0;
        double sum = 0;
        MPI_Reduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        total_sum += sum;
        total_points_c += NUMBER_OF_POINTS;
        free(points);

        result = area * total_sum / total_points_c;
        err = std::abs(result - precise_res);

        if (err <= eps) {
            exiting = 1;
        }
        MPI_Bcast(&exiting, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    std::cout << "Iterations: " << iteration << std::endl;
    std::cout << "Result: " << result << std::endl;
    std::cout << "Error: " << err << std::endl;
}

void worker_routine() {
    int exiting = 0;

    while (!exiting) {
        long bunch_size = NUMBER_OF_POINTS / (n_proc - 1) * 3;
        double *points = (double *) malloc(bunch_size * sizeof(double) + MPI_BSEND_OVERHEAD + 7);

        MPI_Recv(points, bunch_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        double local_sum = 0;
        for (int i = 0; i < bunch_size; i += 3) {
            local_sum += f(points[i], points[i + 1], points[i + 2]);
        }

        MPI_Reduce(&local_sum, nullptr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        MPI_Bcast(&exiting, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Need to pass precision\n");
        exit(1);
    }

    eps = strtod(argv[1], nullptr);

    int err = MPI_Init(&argc, &argv);

    if (err != MPI_SUCCESS)
    {
        fprintf(stderr, "Error while starting! \n");
        MPI_Abort(MPI_COMM_WORLD, err);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

    double start_time = MPI_Wtime();

    if (!rank) {
        master_routine();
    } else {
        worker_routine();
    }

    double end_time = MPI_Wtime();
    double total_time = end_time - start_time;

    double max_total_time;
    MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (!rank) {
        std::cout << "time: " << max_total_time << "s" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
