#define _USE_MATH_DEFINES

#include <cmath>
#include <omp.h>
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

double volume(int dim, int n, int nthreads) {

    int inside = 0;
    #pragma omp parallel num_threads(nthreads) reduction(+:inside)
    {
        int T = omp_get_num_threads();
        int tid = omp_get_thread_num();
        int chunk = n / T;
        if (tid == T-1)
            chunk += n % T;
        int portion = (n / T)*tid;

        unsigned int seed = time(NULL);
        for (int i = 0; i < chunk; i++) {
            //srand(i + portion);
            double s = 0;
            for (int j = 0; j < dim; j++) {
                double x = 2.0 * rand_r(&seed) / RAND_MAX - 1;
                s += x*x;
            }
            if (s <= 1)
                inside += 1;
        }
    }

    double p = 1;
    for (int i = 0; i < dim; i++)
        p *= 2;

    return p * inside / n;

    /*
    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < n; i++) {
        srand(i);
        double s = 0;
        for (int j = 0; j < dim; j++) {
            double x = 2.0 * rand() / RAND_MAX - 1;
            s += x*x;
        }
        if (s <= 1)
            #pragma omp atomic
            inside += 1;
    }
    */
}

int main() {

    cout << "Dim: " << 10 << endl;
    for (int i = 0; i < 10; i++) {
        //double ex = explicit_volume(i+1);

            struct timespec start, finish;
            double elapsed = 0;
            clock_gettime(CLOCK_MONOTONIC, &start);

            double vol = volume(10, 100000000, i+1);

            clock_gettime(CLOCK_MONOTONIC, &finish);
            elapsed = (finish.tv_sec - start.tv_sec);
            elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
            cout << elapsed << " " << i+1 << endl;

        //cout << "Explicit: " << explicit_volume(3) << endl;
    }

    return 0;
}
