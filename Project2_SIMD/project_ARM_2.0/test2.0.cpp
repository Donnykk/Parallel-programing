#include <arm_neon.h>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#define millitime(x) (x.tv_sec * 1000 + x.tv_usec / 1000.0)
using namespace std;
float **arr = NULL;
void initial(int n)
{
	srand(time(0));
	arr = new float *[n];
	for (int i = 0; i < n; i++)
	{
		arr[i] = new float[n];
		for (int j = 0; j < i; j++)
			arr[i][j] = 0;
		arr[i][i] = 1.0;
		for (int j = i + 1; j < n; j++)
			arr[i][j] = rand();
	}
	for (int k = 0; k < n; k++)
		for (int i = k + 1; i < n; i++)
			for (int j = 0; j < n; j++)
				arr[i][j] += arr[k][j];
}

void func(void (*f)(int), int scale)
{
	int counter = 0;
	timeval start, finish, now;
	gettimeofday(&start, NULL);
	gettimeofday(&now, NULL);
	while (millitime(now) - millitime(start) < 10)
	{
		counter++;
		f(scale);
		gettimeofday(&now, NULL);
	}
	gettimeofday(&finish, NULL);
	cout << (millitime(finish) - millitime(start)) / counter << endl;
}

void origin_version(int n)
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

void gauss_1(int n)
{
	for (int k = 0; k < n; k++)
	{
		float index = arr[k][k];
		float32x4_t v1 = vmovq_n_f32(index);
		float32x4_t v0;
		int l;
		for (l = k + 1; l <= n - 4; l += 4)
		{
			v0 = vld1q_f32(arr[k] + l);
			v0 = vdivq_f32(v0, v1);
			vst1q_f32(arr[k] + l, v0);
		}
		for (l; l < n; l++)
			arr[k][l] = arr[k][l] / index;
		arr[k][k] = 1.0;
		for (int i = k + 1; i < n; i++)
		{
			for (int j = k + 1; j < n; j++)
				arr[i][j] = arr[i][j] - arr[i][k] * arr[k][j];
			arr[i][k] = 0;
		}
	}
}

void gauss_2(int n)
{
	for (int k = 0; k < n; k++)
	{
		float index = arr[k][k];
		for (int j = k + 1; j < n; j++)
			arr[k][j] = arr[k][j] / index;
		arr[k][k] = 1.0;
		for (int i = k + 1; i < n; i++)
		{
			float32x4_t v1 = vmovq_n_f32(arr[i][k]);
			float32x4_t v0, v2;
			int j;
			for (j = k + 1; j <= n - 4; j += 4)
			{
				v2 = vld1q_f32(arr[k] + j);
				v0 = vld1q_f32(arr[i] + j);
				v2 = vmulq_f32(v1, v2);
				v0 = vsubq_f32(v0, v2);
				vst1q_f32(arr[i] + j, v0);
			}
			for (j; j < n; j++)
				arr[i][j] = arr[i][j] - arr[i][k] * arr[k][j];
			arr[i][k] = 0;
		}
	}
}

void gauss_plus(int n)
{
	for (int k = 0; k < n; k++)
	{
		float index = arr[k][k];
		float32x4_t v1 = vmovq_n_f32(index);
		float32x4_t v0;
		int l;
		for (l = k + 1; l <= n - 4; l += 4)
		{
			v0 = vld1q_f32(arr[k] + l);
			v0 = vdivq_f32(v0, v1);
			vst1q_f32(arr[k] + l, v0);
		}
		for (l; l < n; l++)
			arr[k][l] = arr[k][l] / index;
		arr[k][k] = 1.0;
		for (int i = k + 1; i < n; i++)
		{
			v1 = vmovq_n_f32(arr[i][k]);
			float32x4_t v2;
			int j;
			for (j = k + 1; j <= n - 4; j += 4)
			{
				v2 = vld1q_f32(arr[k] + j);
				v0 = vld1q_f32(arr[i] + j);
				v2 = vmulq_f32(v1, v2);
				v0 = vsubq_f32(v0, v2);
				vst1q_f32(arr[i] + j, v0);
			}
			for (j; j < n; j++)
				arr[i][j] = arr[i][j] - arr[i][k] * arr[k][j];
			arr[i][k] = 0;
		}
	}
}

int main()
{
	int scale = 64;
	while (scale <= 1024)
	{
		cout << "scale=" << scale << endl;
		initial(scale);
		cout << "origin_version: ";
		func(origin_version, scale);
		cout << "gauss_1: ";
		func(gauss_1, scale);
		cout << "gauss_2: ";
		func(gauss_2, scale);
		cout << "gauss_plus: ";
		func(gauss_plus, scale);
		scale *= 2;
		cout << endl;
	}
	return 0;
}
