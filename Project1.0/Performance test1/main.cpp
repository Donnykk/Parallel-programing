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
        cin >> a[i];
    }
    float** b = new float*[n];
    for(int i=0;i<n;i++){
        b[i] = new float[n];
        for(int j=0;j<n;j++){
            cin >> b[i][j];
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

    float* sum_ = new float[n];
    long long head_, tail_, freq_;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq_);
    QueryPerformanceCounter((LARGE_INTEGER *)&head_);
    for(int i=0;i<n;i++){
        sum_[i] = 0.0;
    }
    for(int j = 0; j < n; j++){
        for(int i = 0; i < n; i++){
            sum_[i] += b[j][i] * a[j];
        }
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail_);
    cout << (tail_ - head_) * 1000.0 / freq_ << "ms" << endl;
    return 0;
}
