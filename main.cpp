#include <iostream>
#include <mpi.h>
#include <cmath>

double a1 = 0, a2 = 0, a3 = 0, b1 = 1, b2 = 1, b3 = 1;

const int NUMBER_OF_POINTS = 1000000;

struct Point {
    double x, y, z;

    Point(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
};

double f(struct Point *point) {
    return sqrt(point->x * point->x + point->y * point->y);
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void master_routine() {
    srand(time(nullptr));

    Point *points = (Point *) malloc(NUMBER_OF_POINTS * sizeof(Point) + MPI_BSEND_OVERHEAD + 7);
    for (int i = 0; i < NUMBER_OF_POINTS; i++) {
        double x = fRand(a1, b1);
        double y = fRand(a2, b2);
        double z = fRand(a3, b3);

        points[i] = Point(x, y, z);
    }

    for (int i = 0; i < 5; i++) {
        std::cout << points[i].x << " " << points[i].y << " " << points[i].z << std::endl;
    }
}

int main(int argc, char *argv[]) {
    int size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (!rank) {
        master_routine();
    } else {
        printf("Hello from %d\n", rank);
        fflush(stdin);
    }

    MPI_Finalize();

    return 0;
}
