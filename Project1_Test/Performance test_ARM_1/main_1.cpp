#include <iostream>
#include<sys/time.h>
#include<math.h>
using namespace std;
void func1(int n)
{
    float* a = new float[n];
    for(int i=0; i<n; i++)
    {
        a[i] = 10 * i;
    }
    float** b = new float*[n];
    for(int i=0; i<n; i++)
    {
        b[i] = new float[n];
        for(int j=0; j<n; j++)
        {
            b[i][j] = i * j;
        }
    }
    float* sum = new float[n];
    struct timeval start,last;
    gettimeofday(&start,NULL);
    for(int t=0; t<1000; t++)
    {
        for(int i=0; i<n; i++)
        {
            sum[i]=0.0;
            for(int j =0; j<n; j++)
            {
                sum [i]+=b[j][i]*a[j];
            }
        }
    }
    gettimeofday(&last,NULL);
    cout << "n = " << n << " :" << endl;
    cout<<"time1= "<<1000000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;


    for(int i=0; i<n; i++)
    {
        sum[i] = 0.0;
    }
    gettimeofday(&start,NULL);
    for(int t=0; t<1000; t++)
    {
        for(int j=0; j<n; j++)
        {
            for(int i=0; i<n; i++)
            {
                sum[i]+=b[j][i]*a[j];
            }
        }
    }
    gettimeofday(&last,NULL);
    cout<<"time2= "<<1000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;
}

void func2(int n)
{
    float* a = new float[n];
    for(int i=0; i<n; i++)
    {
        a[i] = 10 * i;
    }
    float** b = new float*[n];
    for(int i=0; i<n; i++)
    {
        b[i] = new float[n];
        for(int j=0; j<n; j++)
        {
            b[i][j] = i * j;
        }
    }
    float* sum = new float[n];
    struct timeval start,last;
    gettimeofday(&start,NULL);
    for(int i=0; i<n; i++)
    {
        sum[i]=0.0;
        for(int j =0; j<n; j++)
        {
            sum [i]+=b[j][i]*a[j];
        }
    }
    gettimeofday(&last,NULL);
    cout << "n = " << n << " :" << endl;
    cout<<"time1= "<<1000000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;

    gettimeofday(&start,NULL);
    for(int j=0; j<n; j++)
    {
        for(int i=0; i<n; i++)
        {
            sum[i]+=b[j][i]*a[j];
        }
    }
    gettimeofday(&last,NULL);
    cout<<"time2= "<<1000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;
}


int main()
{
    int n=10;
    func1(n);
    n=100;
    func1(n);
    n=1000;
    func2(n);
    n=5000;
    func2(n);
    n=10000;
    func2(n);
    n=15000;
    func2(n);
    n=20000;
    func2(n);
    return 0;
}

