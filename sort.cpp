#include <omp.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fstream>


using std::cout;
using std::endl;
using std::cin;


double make_sorting(long long n, int nthreads, long long* arr) {

    srand(time(NULL));

    long long** parts = new long long*[nthreads];
    long long* what_left = new long long[nthreads];

    struct timespec start, finish;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);

    #pragma omp parallel num_threads(nthreads)
    {
        int T = omp_get_num_threads();
        int tid = omp_get_thread_num();
        long long start = n / T;
        long long chunk = start;
        if (tid == T-1)
            chunk += n % T;

        parts[tid] = new long long[chunk];
        std::copy(arr + tid*start, arr + tid*start + chunk, parts[tid]);

        what_left[tid] = chunk;
        std::sort(parts[tid], parts[tid] + chunk, [](long long x, long long y) { return x > y; });

    }

    for (long long i = 0; i < n; i++) {
        long long min = n*n + 1;
        long long what_j = -1;

        for (long long j = 0; j < nthreads; j++) {
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

    //for (int n = (int)10e7; n < (int)13e7; n += (int)2e7) {
        //for (int nthreads = 1; nthreads <= 8; nthreads++) {
            double* res = new double[5];
            long long N = 10e7;
            long long* arr = new long long[N];

            for (int i = 0; i < 5; i++) {
                std::ifstream fin("arr", std::ifstream::in);
                for (long long i = 0; i < N; i++) {
                    fin >> arr[i];
                }
                cout << endl;
                fin.close();

                double elapsed = make_sorting(N, 5, arr);
                cout << "Sorting took " << elapsed << " secs" << endl;
                res[i] = elapsed;

                bool f = false;
                for (long long i = 0; i < N-1; i++) {
                    //cout << arr[i] << " ";
                    if (arr[i] > arr[i+1]) {
                        cout << "Sorting failed" << endl;
                        f = true;
                        break;
                    }
                }
                if (!f)
                    cout << "Sorting successful" << endl;
                cout << endl;
            }
            delete[] arr;
            for (int i = 0; i < 5; i++)
                cout << res[i] << ",";
            cout << endl;
            delete[] res;
        //}
    //}


    return 0;
}
