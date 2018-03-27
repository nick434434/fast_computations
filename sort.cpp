#include <omp.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <mpi.h>


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
        long long begin = n / T;
        long long chunk = begin;
        if (tid == T-1)
            chunk += n % T;

        parts[tid] = new long long[chunk];
        std::copy(arr + tid*begin, arr + tid*begin + chunk, parts[tid]);

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


void make_sorting(long long n, long long* arr, double* elapsed) {

    MPI_Status status;

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //cout << "MPI initialized" << endl;

    struct timespec start, finish;

    if (rank == 0)
        clock_gettime(CLOCK_MONOTONIC, &start);

    if (size == 1) {
        cout << "Size == 1, no parallel processing" << endl;
        std::sort(arr, arr + n, [](long long x, long long y) { return x < y; });

        clock_gettime(CLOCK_MONOTONIC, &finish);

        *elapsed = (finish.tv_sec - start.tv_sec);
        *elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        return;
    }

    long long** parts = new long long*[size];
    long long* what_left = new long long[size];

    long long begin = n / size;
    long long chunk = begin;
    if (rank == size-1)
        chunk += n % size;

    // parts[rank] = tmp | Send 1
    long long* tmp = new long long[chunk];
    //cout << "Copying " << rank << endl;
    std::copy(arr + rank*begin, arr + rank*begin + chunk, tmp);

    //cout << "Sorting " << rank << endl;
    // what_left[rank] = chunk | Send 2
    std::sort(tmp, tmp + chunk, [](long long x, long long y) { return x > y; });

    if (rank > 0) {
        //cout << "Sending 1" << endl;
        MPI_Send(&chunk, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
        //cout << "Sending 2" << endl;
        MPI_Send(tmp, chunk, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
    } else {

        parts[0] = new long long[chunk];
        parts[0] = tmp;
        what_left[0] = chunk;

        //cout << "Receiving" << endl;
        for (int i = 1; i < size; ++i) {
            //cout << "recv " << i << endl;
            MPI_Recv(what_left + i, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, &status);
            parts[i] = new long long[what_left[i]];
            MPI_Recv(parts[i], what_left[i], MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, &status);
        }

        //cout << "Assigning, what_left:" << endl;
        //for (int i = 0; i < size; i++)
        //    cout << what_left[i] << " ";
        //cout << endl;
        for (long long i = 0; i < n; i++) {
            long long min = n*n + 1;
            long long what_j = -1;

            for (long long j = 0; j < size; j++) {
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

        *elapsed = (finish.tv_sec - start.tv_sec);
        *elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    }

}


int main(int argc, char** argv) {

    MPI_Init(NULL, NULL);

    long long n = 16e7;
    long long* arr = new long long[n];

    unsigned int t = time(NULL);
    for (int j = 0; j < n; j++) {
        arr[j] = (long long)(rand_r(&t) % 1000000);
        // cout << arr[j] << " ";
    }
    cout << endl;

    //cout << "Entering sorting" << endl;
    double elapsed;
    MPI_Status status;

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //cout << "MPI initialized" << endl;

    struct timespec start, finish;

    if (rank == 0)
        clock_gettime(CLOCK_MONOTONIC, &start);

    if (size == 1) {
        cout << "Size == 1, no parallel processing" << endl;
        std::sort(arr, arr + n, [](long long x, long long y) { return x < y; });

        clock_gettime(CLOCK_MONOTONIC, &finish);

        *elapsed = (finish.tv_sec - start.tv_sec);
        *elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        return;
    }

    long long** parts = new long long*[size];
    long long* what_left = new long long[size];

    long long begin = n / size;
    long long chunk = begin;
    if (rank == size-1)
        chunk += n % size;

    // parts[rank] = tmp | Send 1
    long long* tmp = new long long[chunk];
    //cout << "Copying " << rank << endl;
    std::copy(arr + rank*begin, arr + rank*begin + chunk, tmp);

    //cout << "Sorting " << rank << endl;
    // what_left[rank] = chunk | Send 2
    std::sort(tmp, tmp + chunk, [](long long x, long long y) { return x > y; });

    if (rank > 0) {
        //cout << "Sending 1" << endl;
        MPI_Send(&chunk, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
        //cout << "Sending 2" << endl;
        MPI_Send(tmp, chunk, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
    } else {

        parts[0] = new long long[chunk];
        parts[0] = tmp;
        what_left[0] = chunk;

        //cout << "Receiving" << endl;
        for (int i = 1; i < size; ++i) {
            //cout << "recv " << i << endl;
            MPI_Recv(what_left + i, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, &status);
            parts[i] = new long long[what_left[i]];
            MPI_Recv(parts[i], what_left[i], MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, &status);
        }

        //cout << "Assigning, what_left:" << endl;
        //for (int i = 0; i < size; i++)
        //    cout << what_left[i] << " ";
        //cout << endl;
        for (long long i = 0; i < n; i++) {
            long long min = n*n + 1;
            long long what_j = -1;

            for (long long j = 0; j < size; j++) {
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

        *elapsed = (finish.tv_sec - start.tv_sec);
        *elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

        cout << elapsed << endl;

        bool f = false;
        for (long long i = 0; i < N-1; i++) {
            // cout << arr[i] << " ";
            if (arr[i] > arr[i+1]) {
                cout << "Sorting failed" << endl;
                f = true;
                break;
            }
        }
        // cout << arr[N-1] << endl;
        if (!f)
            cout << "Sorting successful" << endl;
        cout << endl;
        delete[] arr;

    }

    MPI_Finalize();

    return 0;
}
