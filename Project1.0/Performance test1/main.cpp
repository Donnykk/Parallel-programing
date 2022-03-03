#include <iostream>
#include <windows.h>
#include <stdlib.h>
using namespace std;

int main()
{
    int matrix[5][5] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
    int vec[5] = {10,20,30,40,50};
    int result[5] = {0};
    long long head, tail, freq;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    //ÖðÁÐ¼ÆËã
    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++){
            result[i] += matrix[j][i] * vec[j];
        }
        cout << "col " << i+1 << ":" << result[i] << endl;
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    cout << (tail - head) * 1000.0 / freq << "ms" << endl;
    return 0;
}
