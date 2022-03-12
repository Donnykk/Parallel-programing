#include <iostream>
#include<sys/time.h>
#include<math.h>
using namespace std;
void func(double* a)
{
    for(int m=10000; m>1; m/=2)
    {
        for(int i=0; i<m/2; i++)
        {
            a[i]=a[i*2]+a[i*2+1];
        }
    }
}
int main()
{
    int n = 0;
    cin >> n;
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
    cout<<"time= "<<1000000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;

    int sum1=0,sum2=0;
    gettimeofday(&start,NULL);
    for(int i=0;i<1000;i++)
    {
        sum=0;
        sum1=0;
        sum2=0;
        for(int i=0;i<n;i+=2){
            sum1+=a[i];
            sum2+=a[i+1];
        }
        sum = sum1+sum2;
    }
    gettimeofday(&last,NULL);
    cout<<"time= "<<1000000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;

    gettimeofday(&start,NULL);
    for(int i=0; i<1000; i++)
    {
        func(a);
    }
    gettimeofday(&last,NULL);
    cout<<"time= "<<1000000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;
    return 0;
}
