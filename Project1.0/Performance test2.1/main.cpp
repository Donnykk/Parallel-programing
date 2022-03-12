#include <iostream>
#include <windows.h>
using namespace std;
int main()
{
    int n = 0, sum = 0;
    cin >> n;
    int* arr = new int[n];
    for(int i=0;i<n;i++){
        arr[i] = i+100;
    }

    int sum1 = 0, sum2 = 0;
    long long head, tail, freq;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int i=0;i<n;i+=2){
        sum1+=arr[i];
        sum2+=arr[i+1];
    }
    sum = sum1+sum2;
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    cout << "sum = " << sum << "  " << (tail - head) * 1000.0 / freq << "ms" << endl;
    return 0;
}
