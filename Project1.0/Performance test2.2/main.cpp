#include <iostream>
#include <windows.h>
using namespace std;

int sum_(int* arr, int n){
    if(n==1){
        return arr[0];
    } else {
        for(int i=0;i<n/2;i++){
            arr[i] += arr[n-i-1];
        }
        n = n/2;
        sum_(arr, n);
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
    sum = sum_(arr, n);
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    cout << "sum = " << sum << "  " << (tail - head) * 1000.0 / freq << "ms" << endl;
    return 0;
}
