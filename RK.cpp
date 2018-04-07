#define _USE_MATH_DEFINES

#include <cmath>
#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>


using std::cout;
using std::endl;


double rhs(double* xy, double* rxy, double* omega, int n, int k, int nthreads) {
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

def RK(int niter, int n, int nthreads, int k, double h, double* xy0, double** res, double* w) {

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
        memset(res[i], 0, nn*4*sizeof(double));

        double* last = res[i-1];

        rhs(res[i-1], k1, w, n, k, nthreads);

        for (int j = 0; j < n; ++j)
            tmp[j] = last[j] + k1[j] * h_2;
        rhs(tmp, k2, w, n, k, nthreads);

        for (int j = 0; j < n; ++j)
            tmp[j] = last[j] + k2[j] * h_2;
        rhs(tmp, k3, w, n, k, nthreads);

        for (int j = 0; j < n; ++j)
            tmp[j] = last[j] + k3[j] * h_2;
        rhs(tmp, k4, w, n, k, nthreads);

        for (int j = 0; j < n; ++j)
            res[i][j] = last[j] + h_6 * (k1[j] + k4[j] + 2 * (k2[j] + k3[j]));

    }

    delete[] k1;
    delete[] tmp;
}


int main(int argc, char** argv) {



    return 0;
}
