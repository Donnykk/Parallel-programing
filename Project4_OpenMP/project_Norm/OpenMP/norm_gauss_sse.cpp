#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <immintrin.h>
#include <sys/time.h>
#define millitime(x) (x.tv_sec * 1000 + x.tv_usec / 1000.0)
using namespace std;
float **arr = NULL;
int count_T = 8, n = 64;

void initial(int n)
{
	float **temp = new float *[n];
	arr = new float *[n];
	for (int i = 0; i < n; i++)
	{
		arr[i] = new float[n]{0};
		temp[i] = new float[n]{0};
	}
	//避免inf和nan,先初始化为上三角矩阵
	for (int i = 0; i < n; i++)
	{
		int j = i;
		arr[i][j] = 1;
		temp[i][j++] = 1;
		for (; j < n; j++)
		{
			arr[i][j] = i + j;
			temp[i][j] = arr[i][j];
		}
	}
	for (int i = 1; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				arr[i][k] += temp[j][k];
}

void origin_version()
{
	for (int k = 0; k < n; k++)
	{
		for (int j = k + 1; j < n; j++)
			arr[k][j] /= arr[k][k];
		arr[k][k] = 1.0;
		for (int i = k + 1; i < n; i++)
		{
			for (int j = k + 1; j < n; j++)
				arr[i][j] = arr[i][j] - arr[i][k] * arr[k][j];
			arr[i][k] = 0;
		}
	}
}

void gauss_OpenMP()
{
	int i, j, k;
	float temp;
#pragma omp parallel num_threads(count_T), private(i, j, k, temp)
	for (k = 0; k < n; k++)
	{
//串行部分
#pragma omp single
		{
			temp = arr[k][k];
			for (j = k + 1; j < n; j++)
				arr[k][j] = arr[k][j] / temp;
			arr[k][k] = 1.0;
		}
//并行部分
#pragma omp for
		for (i = k + 1; i < n; i++)
		{
			temp = arr[i][k];
			for (j = k + 1; j < n; j++)
				arr[i][j] = arr[i][j] - temp * arr[k][j];
			arr[i][k] = 0;
		}
	}
}

void gauss_SSE()
{
	__m128 vt, va, temp1, temp2, temp3;
	int i, j, k;
	float temp;
#pragma omp parallel num_threads(count_T), private(i, j, k, temp)
	for (k = 0; k < n; k++)
	{
#pragma omp single
		{
			temp = arr[k][k];
			for (j = k + 1; j < n; j++)
				arr[k][j] = arr[k][j] / temp;
			arr[k][k] = 1.0;
		}
#pragma omp for
		for (i = k + 1; i < n; i++)
		{
			for (j = k; j < n; j += 4)
			{
				if (j + 4 > n)
				{
					for (; j < n; j++)
						arr[i][j] = arr[i][j] - arr[i][k] * arr[k][j];
				}
				else
				{
					temp1 = _mm_loadu_ps(arr[i] + j);
					temp2 = _mm_loadu_ps(arr[k] + j);
					temp3 = _mm_set1_ps(arr[i][k]);
					temp2 = _mm_mul_ps(temp3, temp2);
					temp1 = _mm_sub_ps(temp1, temp2);
					_mm_storeu_ps(arr[i] + j, temp1);
				}
				arr[i][k] = 0;
			}
		}
	}
}

void func(void (*f)())
{
	int counter = 0;
	initial(n);
	timeval start, finish, now;
	gettimeofday(&start, NULL);
	gettimeofday(&now, NULL);
	while (millitime(now) - millitime(start) < 10)
	{
		counter++;
		f();
		gettimeofday(&now, NULL);
	}
	gettimeofday(&finish, NULL);
	cout << (millitime(finish) - millitime(start)) / counter << endl;
}

int main()
{
	while (n <= 2048)
	{
		cout << "scale=" << n << endl;
		//平凡算法
		cout << endl
			 << "origin_version: ";
		func(origin_version);
		// OpenMP优化
		cout << endl
			 << "gauss_OpenMP: ";
		func(gauss_OpenMP);
		// 结合SSE
		cout << endl
			 << "gauss_SSE: ";
		func(gauss_SSE);
		n *= 2;
		cout << endl;
	}
	return 0;
}
