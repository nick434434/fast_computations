#include <omp.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>


using std::cout;
using std::endl;
using std::cin;


double make_sorting(int n, int nthreads, int* arr) {

    srand(time(NULL));

    cout << n << " " << nthreads << endl;

    for (int i = 0; i < n; i++) {
        int a = rand();
        //cout << a << ":";
        arr[i] = (a / n);
        //cout << arr[i] << " ";
    }

    int** parts = new int*[nthreads];
    int* what_left = new int[nthreads];

    struct timespec start, finish;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);

    #pragma omp parallel num_threads(nthreads)
    {
        int T = omp_get_num_threads();
        int tid = omp_get_thread_num();
        int start = n / T;
        int chunk = start;
        if (tid == T-1)
            chunk += n % T;

        parts[tid] = new int[chunk];
        std::copy(arr + tid*start, arr + tid*start + chunk, parts[tid]);

        what_left[tid] = chunk;
        std::sort(parts[tid], parts[tid] + chunk, [](int x, int y) { return x > y; });

    }

    for (int i = 0; i < n; i++) {
        int min = n*n + 1;
        int what_j = -1;

        for (int j = 0; j < nthreads; j++) {
            if (what_left[j]) {
                if (parts[j][what_left[j]-1] < min) {
                    min = parts[j][what_left[j]-1];
                    what_j = j;
                }
            }
        }
        what_left[what_j] -= 1;
        arr[i] = min;
    }

    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    return elapsed;
}


int main(int argc, char** argv) {

    for (int n = 10000000; n < 200000000; n *= 2) {
        for (int nthreads = 1; nthreads <= 8; nthreads++) {
            int* arr = new int[n];
            double elapsed = make_sorting(n, nthreads, arr);
            cout << "Sorting took " << elapsed << " secs" << endl;

            bool f = false;
            for (int i = 0; i < n-1; i++) {
                //cout << arr[i] << " ";
                if (arr[i] > arr[i+1]) {
                    cout << "Sorting failed" << endl;
                    f = true;
                    break;
                }
            }
            delete[] arr;
            if (!f)
                cout << "Sorting successful" << endl;
            cout << endl;
        }
    }


    return 0;
}
