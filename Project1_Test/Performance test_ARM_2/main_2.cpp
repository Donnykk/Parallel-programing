#include <iostream>
#include<sys/time.h>
#include<math.h>
using namespace std;
void recursion(double *a,int n)
{
    if(n==1)
    {
        return;
    }
    else
    {
        for(int i=0;i<n/2;i++)
        {
            a[i] += a[n-i-1];
        }
        n/=2;
        recursion(a,n);
    }
}

void func1(int n)
{
    double* a = new double[n];
    for(int i=0; i<n; i++)
    {
        a[i] = 10*i;
    }
    double sum = 0;
    struct timeval start,last;
    gettimeofday(&start,NULL);
    for(int cnt=0; cnt<1000; cnt++)
    {
        sum=0;
        for(int i=0; i<n; i++)
            sum+=a[i];
    }
    gettimeofday(&last,NULL);
    cout << "n = " << n << " :" << endl;
    cout<<"time1= "<<1000000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;

    int sum1=0,sum2=0;
    gettimeofday(&start,NULL);
    for(int i=0; i<1000; i++)
    {
        sum=0;
        sum1=0;
        sum2=0;
        for(int i=0; i<n; i+=2)
        {
            sum1+=a[i];
            sum2+=a[i+1];
        }
        sum = sum1+sum2;
    }
    gettimeofday(&last,NULL);
    cout<<"time2= "<<1000000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;

    gettimeofday(&start,NULL);
    for(int i=0; i<1000; i++)
    {
        double* b = new double[n];
        for(int i=0;i<n;i++)
            b[i] = a[i];
        recursion(b,n);
    }
    gettimeofday(&last,NULL);
    cout<<"time3= "<<1000000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;
}

void func2(int n)
{
    double* a = new double[n];
    for(int i=0; i<n; i++)
    {
        a[i] = 10*i;
    }
    double sum = 0;
    struct timeval start,last;
    gettimeofday(&start,NULL);
    for(int i=0; i<n; i++)
        sum+=a[i];
    gettimeofday(&last,NULL);
    cout << "n = " << n << " :" << endl;
    cout<<"time1= "<<1000000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;

    int sum1=0,sum2=0;
    gettimeofday(&start,NULL);
    sum=0;
    for(int i=0; i<n; i+=2)
    {
        sum1+=a[i];
        sum2+=a[i+1];
    }
    sum = sum1+sum2;
    gettimeofday(&last,NULL);
    cout<<"time2= "<<1000000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;

    gettimeofday(&start,NULL);
    double* b = new double[n];
    for(int i=0;i<n;i++)
        b[i] = a[i];
    recursion(b,n);
    gettimeofday(&last,NULL);
    cout<<"time3= "<<1000000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;
}

int main()
{
    int n = 256;
    func1(n);
    n = 1024;
    func1(n);
    n = 8192;
    func1(n);
    n = 16384;
    func1(n);
    n = 65536;
    func1(n);
    n = 131072;
    func2(n);
    n = 262144;
    func2(n);
    return 0;
}
