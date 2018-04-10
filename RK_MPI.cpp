#define _USE_MATH_DEFINES

#include <cmath>
#include <omp.h>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <cstring>
#include <mpi.h>
#include <fstream>


using std::cout;
using std::endl;


void rhs(double* xy, double* rxy, double* omega, int n, int k, int nthreads) {
    double* x = xy;
    double* y = xy + n;
    double* rx = rxy;
    double* ry = rxy + n;

    int cluster_size = n / nthreads;
    int begin = 0, end = 0;
    if (k < n % nthreads) {
        cluster_size += 1;
        begin += k;
        end += k;
    } else {
        begin  += n % nthreads;
        end += n % nthreads;
    }

    begin += k * cluster_size;
    end += (k + 1) * cluster_size;

    for (int i = 0; i < n; i++) {
        for (int j = begin; j < end; j++) {
            if (i == j)
                continue;

            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double t = omega[j] / (dx*dx + dy*dy);
            rx[i] += dy * t;
            ry[i] += dx * t;
        }
        rx[i] /= 2*M_PI;
        ry[i] /= 2*M_PI;
    }

    MPI_Allreduce(MPI_IN_PLACE, rxy, 2*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void RK(int niter, int n, int nthreads, int k, double h, double* xy0, double**& res, double* w) {

    int nn = n*2;
    double h_2 = h / 2.;
    double h_6 = h / 6.;

    double* k1 = new double[nn*4];
    double* k2 = k1 + nn;
    double* k3 = k2 + nn;
    double* k4 = k3 + nn;
    double* tmp = new double[nn];

    res = new double*[niter];

    res[0] = new double[nn];
    memcpy(res[0], xy0, nn*sizeof(double));

    for (int i = 1; i < niter; ++i) {
        memset(k1, 0, nn*4*sizeof(double));
        res[i] = new double[nn];
        memset(res[i], 0, nn*sizeof(double));

        double* last = res[i-1];
        rhs(res[i-1], k1, w, n, k, nthreads);

        for (int j = 0; j < nn; ++j)
            tmp[j] = last[j] + k1[j] * h_2;
        rhs(tmp, k2, w, n, k, nthreads);

        for (int j = 0; j < nn; ++j)
            tmp[j] = last[j] + k2[j] * h_2;
        rhs(tmp, k3, w, n, k, nthreads);

        for (int j = 0; j < nn; ++j)
            tmp[j] = last[j] + k3[j] * h_2;
        rhs(tmp, k4, w, n, k, nthreads);

        for (int j = 0; j < nn; ++j)
            res[i][j] = last[j] + h_6 * (k1[j] + k4[j] + 2 * (k2[j] + k3[j]));

    }

    delete[] k1;
    delete[] tmp;
}


void gen(double* xy, double* w, int n) {

    unsigned int seed = (unsigned int)time(NULL);
    for (int i = 0; i < n; i++) {
        double arg = 1.0 * rand_r(&seed) / RAND_MAX;
        arg *= 2.0 * M_PI;
        xy[i] = cos(arg);
        xy[i + n] = sin(arg);
        w[i] = 1;
    }

}

void save(std::string fname, double time, double** res, double* w, int niter, int n, double h) {

    std::ofstream fout(fname);

    fout << time << "\n";
    fout << n << " " << niter << " " << h << "\n";

    for (int i = 0; i < n; i++)
        fout << w[i] << " ";
    fout << "\n";

    for (int i = 0; i < niter; i++) {
        double* r = res[i];
        for (int j = 0; j < 2*n; j++) {
            if (j == n)
                fout << "\n";
            fout << r[j] << " ";
        }
        fout << "\n";

    }

}

void run(int n, int niter) {

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double* xy = new double[2*n];
    double* w = new double[n];
    if (rank == 0) {
        gen(xy, w, n);
    }

    MPI_Bcast(xy, 2*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(w, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < 2*n; i++)
            cout << xy[i] << " ";
        cout << endl;
    }

    double** res;
    double h = 0.01;

    struct timespec start, finish;
    double elapsed = 0;
    if (rank == 0) {
        clock_gettime(CLOCK_MONOTONIC, &start);
    }

    RK(niter, n, size, rank, h, xy, res, w);

    if (rank == 0) {
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        save("points", elapsed, res, w, niter, n, h);
    }

}


int main(int argc, char** argv) {

    MPI_Init(NULL, NULL);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
        cout << argc << endl;
    run(atoi(argv[1]), atoi(argv[2]));

    MPI_Finalize();

    return 0;
}
