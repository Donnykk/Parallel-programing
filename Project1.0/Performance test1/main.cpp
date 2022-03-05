#include <iostream>
#include <windows.h>
#include <stdlib.h>
using namespace std;

int main()
{
    int n = 0;
    cin >> n;
    float* a = new float[n];
    for(int i=0;i<n;i++){
        a[i] = 10 * i;
    }
    float** b = new float*[n];
    for(int i=0;i<n;i++){
        b[i] = new float[n];
        for(int j=0;j<n;j++){
            b[i][j] = i * j;
        }
    }
    float* sum = new float[n];
    long long head, tail, freq;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int i=0;i<n;i++){
        sum[i] = 0.0;
        for(int j=0;j<n;j++){
            sum[i] += b[j][i] * a[j];
        }
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    cout << (tail - head) * 1000.0 / freq << "ms" << endl;
    return 0;
}
