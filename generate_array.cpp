#include <iostream>
#include <fstream>
#include <ctime>

using std::cout;
using std::endl;

int main() {

    const long long N = 10e7;
    const uint k = 5;

    std::ofstream fout("arr", std::ofstream::out | std::ofstream::trunc);
    //int* a = new int[N];
    unsigned int seed = time(NULL);
    for (long long i = 0; i < N; i++) {
        //a[i] = i;
        fout << rand_r(&seed) << " ";
    }
    cout << endl;

    fout.close();

    return 0;
}
