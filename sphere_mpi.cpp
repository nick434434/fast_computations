#define _USE_MATH_DEFINES

#include <cmath>
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>


using std::cout;
using std::endl;
using std::cin;


double explicit_volume(int dim) {
    if (dim == 1)
        return 2;
    if (dim == 2)
        return M_PI;
    return 2 * M_PI / dim * explicit_volume(dim - 2);
}

int main(int argc, char** argv) {

    MPI_Init(NULL, NULL);
    MPI_Status status;

    struct timespec start, finish;
    double elapsed = 0;
    clock_gettime(CLOCK_MONOTONIC, &start);

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        cout << argc << endl;

    int n = atoi(argv[1]);
    int dim = atoi(argv[2]);
    int inside = 0;
    int chunk = n / size;
    if (rank == size-1)
        chunk += n % size;

    unsigned int seed = time(NULL);
    double s = 0, x = 0;
    for (int i = 0; i < chunk; i++) {
        //srand(i + portion);
        s = 0;
        for (int j = 0; j < dim; j++) {
            x = 2.0 * rand_r(&seed) / RAND_MAX - 1;
            s += x*x;
        }
        if (s <= 1)
            inside += 1;
    }

    if (rank > 0)
        MPI_Send(&inside, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    else {

        int insides;
        for (int i = 1; i < size; i++) {
            MPI_Recv(&insides, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            inside += insides;
        }

        double p = 1;
        for (int i = 0; i < dim; i++)
            p *= 2;

        cout << p * inside / n << " " << dim << " " << n << endl;
    }

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    cout << elapsed << endl;

    MPI_Finalize();

    return 0;
}
