#include <arm_neon.h>
#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#define millitime(x) (x.tv_sec * 1000 + x.tv_usec / 1000.0)
using namespace std;
float **arr = NULL;
int count_T = 4, n = 64;

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
//并行部分c
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

void gauss_Neon()
{
	int i, j, k;
	int chunksize = 10;
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
			float32x4_t v1 = vmovq_n_f32(arr[i][k]);
			float32x4_t v0, v2;
			for (j = k + 1; j <= n - 4; j += 4)
			{
				v2 = vld1q_f32(arr[k] + j);
				v0 = vld1q_f32(arr[i] + j);
				v2 = vmulq_f32(v1, v2);
				v0 = vsubq_f32(v0, v2);
				vst1q_f32(arr[i] + j, v0);
			}
			for (; j < n; j++)
				arr[i][j] = arr[i][j] - arr[i][k] * arr[k][j];
			arr[i][k] = 0;
		}
	}
}

void guass_OpenMp_staticChunk()
{
	int i, j, k;
	int chunksize = 128;
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
#pragma omp for schedule(static, chunksize)
		for (i = k + 1; i < n; i++)
		{
			temp = arr[i][k];
			for (j = k + 1; j < n; j++)
				arr[i][j] = arr[i][j] - temp * arr[k][j];
			arr[i][k] = 0;
		}
	}
}

void guass_OpenMp_dynamicChunk()
{
	int i, j, k;
	int chunksize = 128;
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
#pragma omp for schedule(dynamic, chunksize)
		for (i = k + 1; i < n; i++)
		{
			temp = arr[i][k];
			for (j = k + 1; j < n; j++)
				arr[i][j] = arr[i][j] - temp * arr[k][j];
			arr[i][k] = 0;
		}
	}
}

void guass_OpenMp_guidedChunk()
{
	int i, j, k;
	int chunksize = 128;
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
#pragma omp for schedule(guided, chunksize)
		for (i = k + 1; i < n; i++)
		{
			temp = arr[i][k];
			for (j = k + 1; j < n; j++)
				arr[i][j] = arr[i][j] - temp * arr[k][j];
			arr[i][k] = 0;
		}
	}
}

void func(void (*f)())
{
	int counter = 0;
	timeval begin, start, finish, now;
	float milliseconds = 0;
	gettimeofday(&begin, NULL);
	gettimeofday(&now, NULL);
	while (millitime(now) - millitime(begin) < 10)
	{
		initial(n);
		counter++;
		gettimeofday(&start, NULL);
		f();
		gettimeofday(&finish, NULL);
		milliseconds += (millitime(finish) - millitime(start));
		gettimeofday(&now, NULL);
	}
	cout << milliseconds / counter << endl;
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
		// 结合Neon
		cout << endl
			 << "gauss_Neon: ";
		func(gauss_Neon);
		/*
		// static线程划分
		cout << endl
			 << "gauss_OpenMP_staticChunk: ";
		func(guass_OpenMp_staticChunk);
		// dynamic线程划分
		cout << endl
			 << "gauss_OpenMP_dynamicChunk: ";
		func(guass_OpenMp_dynamicChunk);
		// guided线程划分
		cout << endl
			 << "gauss_OpenMP_guidedChunk: ";
		func(guass_OpenMp_guidedChunk);
		*/
		n *= 2;
		cout << endl;
	}
	return 0;
}
