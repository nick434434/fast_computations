#include <omp.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>


using std::cout;
using std::endl;
using std::cin;


int main(int argc, char** argv) {

    if (argc != 3)
        return -1;

    int nthreads;

    srand(time(NULL));

    int n = atoi(argv[1]);
    nthreads = atoi(argv[2]);
    cout << n << endl;
    cout << nthreads << endl;

    int* arr = new int[n];
    for (int i = 0; i < n; i++) {
        int a = rand();
        cout << a << ":";
        arr[i] = (a % n*n) + 1;
        cout << arr[i] << " ";
    }
    cout << endl;

    int** parts = new int*[nthreads];
    int* what_left = new int[nthreads];

    #pragma omp parallel num_threads(nthreads)
    {
        int T = omp_get_num_threads();
        int tid = omp_get_thread_num();
        int start = n / T;
        int chunk = start;
        if (tid == T-1)
            chunk += n % T;

        #pragma omp critical
        {
            cout << tid*start << " " << chunk << endl;
        }

        parts[tid] = new int[chunk];
        std::copy(arr + tid*start, arr + tid*start + chunk, parts[tid]);

        what_left[tid] = chunk;
        std::sort(parts[tid], parts[tid] + chunk, [](int x, int y) { return x > y; });
        #pragma omp critical
        {
            for (int j = 0; j < chunk; j++) {
                cout << parts[tid][j] << " ";
            }
            cout << endl;
        }
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

    for (int i = 0; i < n-1; i++) {
        cout << arr[i] << " ";
        if (arr[i] > arr[i+1]) {
            cout << "Sorting failed" << endl;
            return 0;
        }
    }

    cout << "Sorting successful" << endl;
    return 0;
}
