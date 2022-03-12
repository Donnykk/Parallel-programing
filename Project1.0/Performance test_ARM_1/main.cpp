#include <iostream>
#include<sys/time.h>
#include<math.h>
using namespace std;
int main()
{
    int n=0;
    cin>>n;
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
    cout<<"time= "<<1000000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;

    gettimeofday(&start,NULL);
    for(int j=0; j<n; j++)
    {
        for(int i=0; i<n; i++)
        {
            sum[i]+=b[j][i]*a[j];
        }
    }
    gettimeofday(&last,NULL);
    cout<<"time= "<<1000*(last.tv_sec-start.tv_sec)+(last.tv_usec-start.tv_usec) << "ms" <<endl;
    return 0;
}

