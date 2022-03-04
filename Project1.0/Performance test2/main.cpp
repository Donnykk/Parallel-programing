#include <iostream>
#include <windows.h>
using namespace std;

int sum_(int* arr, int n){
    if(n==1){
        return arr[0];
    } else {
        if(n % 2 == 0){
            for(int i=0;i<n/2;i++){
                arr[i] += arr[n-i-1];
            }
            n = n/2;
            sum_(arr, n);
        } else {
            for(int i=0;i<n/2;i++){
                arr[i] += arr[n-i-1];
            }
            n = n/2 + 1;
            sum_(arr, n);
        }
    }
}

int main()
{
    int n = 0, sum = 0;
    cin >> n;
    int* arr = new int[n];
    for(int i=0;i<n;i++){
        arr[i] = 10*i + 10;
    }

    long long head, tail, freq;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int i=0;i<n;i++){
        sum+=arr[i];
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    cout << "sum1 = " << sum << "  " << (tail - head) * 1000.0 / freq << "ms" << endl;

    int sum1 = 0, sum2 = 0;
    long long head_, tail_, freq_;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq_);
    QueryPerformanceCounter((LARGE_INTEGER *)&head_);
    if(n % 2 == 0){
        for(int i=0;i<n;i+=2){
            sum1+=arr[i];
            sum2+=arr[i+1];
        }
    } else{
         for(int i=0;i<n-1;i+=2){
            sum1+=arr[i];
            sum2+=arr[i+1];
        }
        sum1+=arr[n-1];
    }
    sum = sum1+sum2;
    QueryPerformanceCounter((LARGE_INTEGER *)&tail_);
    cout << "sum2 = " << sum << "  " << (tail_ - head_) * 1000.0 / freq_ << "ms" << endl;

    long long head_new, tail_new, freq_new;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq_new);
    QueryPerformanceCounter((LARGE_INTEGER *)&head_new);
    sum = sum_(arr, n);
    QueryPerformanceCounter((LARGE_INTEGER *)&tail_new);
    cout << "sum3 = " << sum << "  " << (tail_new - head_new) * 1000.0 / freq_new << "ms" << endl;
    return 0;
}
