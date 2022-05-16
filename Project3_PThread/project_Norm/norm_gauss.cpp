#include <arm_neon.h>
#include <semaphore.h>
#include <pthread.h>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#define millitime(x) (x.tv_sec * 1000 + x.tv_usec / 1000.0)
using namespace std;
float **arr = NULL;
int count_T = 8, n = 64;

typedef struct
{
	int k;	  //消去的轮次
	int t_id; //线程id
} threadParam_td;

typedef struct
{
	int t_id;
} threadParam_ts;

pthread_barrier_t barrier_Divsion;
pthread_barrier_t barrier_Elimination;

sem_t *sem_Division;
sem_t *sem_Elimination;
sem_t sem_leader;

void *threadFunc_dynamic(void *param)
{
	threadParam_td *p = (threadParam_td *)param;
	int k = p->k;
	int i = k + p->t_id + 1;
	for (; i < n; i += count_T)
	{
		for (int j = k + 1; j < n; j++)
			arr[i][j] = arr[i][j] - arr[i][k] * arr[k][j];
		arr[i][k] = 0;
	}
	pthread_exit(NULL);
}

void *threadFunc_static(void *param)
{
	long long t_id = (long long)param;
	for (int k = 0; k < n; k++)
	{
		if (t_id == 0)
		{
			for (int j = k + 1; j < n; j++)
				arr[k][j] = arr[k][j] / arr[k][k];
			arr[k][k] = 1.0;
		}
		pthread_barrier_wait(&barrier_Divsion); //第一个同步点
		//循环划分任务
		for (int i = k + 1 + t_id; i < n; i += count_T)
		{
			for (int j = k + 1; j < n; j++)
				arr[i][j] = arr[i][j] - arr[i][k] * arr[k][j];
			arr[i][k] = 0.0;
		}
		pthread_barrier_wait(&barrier_Elimination);
	}
	pthread_exit(NULL);
}

//结合SIMD
void *threadFunc_SIMD_static(void *param)
{
	long long t_id = (long long)param;
	for (int k = 0; k < n; k++)
	{
		if (t_id == 0)
		{
			for (int j = k + 1; j < n; j++)
				arr[k][j] = arr[k][j] / arr[k][k];
			arr[k][k] = 1.0;
		}
		pthread_barrier_wait(&barrier_Divsion);
		for (int i = k + 1 + t_id; i < n; i += count_T)
		{
			float32x4_t v1 = vmovq_n_f32(arr[i][k]);
			float32x4_t v0, v2;
			int j = k + 1;
			for (; j <= n - 4; j += 4)
			{
				v2 = vld1q_f32(arr[k] + j);
				v0 = vld1q_f32(arr[i] + j);
				v2 = vmulq_f32(v1, v2);
				v0 = vsubq_f32(v0, v2);
				vst1q_f32(arr[i] + j, v0);
			}
			for (; j < n; j++)
				arr[i][j] = arr[i][j] - arr[i][k] * arr[k][j];
			arr[i][k] = 0.0;
		}
		pthread_barrier_wait(&barrier_Elimination);
	}
	pthread_exit(NULL);
	return NULL;
}

void *threadFunc_col(void *param)
{
	long long t_id = (long long)param;
	for (int k = 0; k < n; k++)
	{
		if (t_id == 0)
		{
			for (int j = k + 1; j < n; j++)
				arr[k][j] = arr[k][j] / arr[j][j];
			arr[k][k] = 1.0;
			for (int i = 0; i < count_T - 1; i++)
				sem_post(sem_Division + i);
		}
		else
			sem_wait(sem_Division + t_id - 1);
		int interval = ((n - k - 1) % count_T == 0) ? (n - k - 1) / count_T : (n - k - 1) / count_T + 1;
		int upp_bound = min((long long)n, k + 1 + interval * (t_id + 1));
		for (int i = k + 1; i < n; i++)
		{
			for (int j = k + t_id * interval + 1; j < upp_bound; j++)
				arr[i][j] = arr[i][j] - arr[i][k] * arr[k][j];
		}
		if (t_id == 0)
		{
			for (int i = 0; i < count_T - 1; i++)
				sem_wait(&sem_leader);
			for (int i = k + 1; i < n; i++)
				arr[i][k] = 0;
			for (int i = 0; i < count_T - 1; i++)
				sem_post(sem_Elimination + i);
		}
		else
		{
			sem_post(&sem_leader);
			sem_wait(sem_Elimination + t_id - 1);
		}
	}
	pthread_exit(NULL);
	return NULL;
}

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

//动态线程版本
void gauss_dynamic()
{
	for (int k = 0; k < n; k++)
	{
		for (int j = k + 1; j < n; j++)
			arr[k][j] = arr[k][j] / arr[k][k];
		arr[k][k] = 1.0;
		int worker_count = count_T;
		pthread_t *handles = new pthread_t[worker_count];
		threadParam_td *param = new threadParam_td[worker_count];
		for (int t_id = 0; t_id < worker_count; t_id++)
		{
			param[t_id].k = k;
			param[t_id].t_id = t_id;
		}
		for (int t_id = 0; t_id < worker_count; t_id++)
			pthread_create(handles + t_id, NULL, threadFunc_dynamic, param + t_id);
		for (int t_id = 0; t_id < worker_count; t_id++)
			pthread_join(handles[t_id], NULL);
		delete[] handles;
		delete[] param;
	}
}

//静态线程版本
void gauss_static()
{
	pthread_barrier_init(&barrier_Divsion, NULL, count_T);
	pthread_barrier_init(&barrier_Elimination, NULL, count_T);
	//创建线程
	pthread_t *handles = new pthread_t[count_T];
	threadParam_ts *param = new threadParam_ts[count_T];
	for (int t_id = 0; t_id < count_T; t_id++)
	{
		param[t_id].t_id = t_id;
		pthread_create(handles + t_id, NULL, threadFunc_static, (void *)(long long)t_id);
	}
	for (int i = 0; i < count_T; i++)
		pthread_join(handles[i], NULL);
	delete[] handles;
	delete[] param;
	pthread_barrier_destroy(&barrier_Divsion);
	pthread_barrier_destroy(&barrier_Elimination);
}

void gauss_SIMD_static()
{
	pthread_barrier_init(&barrier_Divsion, NULL, count_T);
	pthread_barrier_init(&barrier_Elimination, NULL, count_T);
	pthread_t *handles = new pthread_t[count_T];
	threadParam_ts *param = new threadParam_ts[count_T];
	for (int t_id = 0; t_id < count_T; t_id++)
	{
		param[t_id].t_id = t_id;
		pthread_create(handles + t_id, NULL, threadFunc_SIMD_static, (void *)(long long)t_id);
	}
	for (int i = 0; i < count_T; i++)
		pthread_join(handles[i], NULL);
	delete[] handles;
	delete[] param;
	pthread_barrier_destroy(&barrier_Divsion);
	pthread_barrier_destroy(&barrier_Elimination);
}

void gauss_col()
{
	sem_Division = new sem_t[count_T - 1];
	sem_Elimination = new sem_t[count_T - 1];
	sem_init(&sem_leader, 0, 0);
	for (int i = 0; i < count_T - 1; i++)
	{
		sem_init(sem_Division + i, 0, 0);
		sem_init(sem_Elimination + i, 0, 0);
	}
	pthread_t *handles = new pthread_t[count_T];
	for (int i = 0; i < count_T; i++)
		pthread_create(handles + i, NULL, threadFunc_col, (void *)(long long)i);
	for (int i = 0; i < count_T; i++)
		pthread_join(handles[i], NULL);
	delete[] handles;
	for (int i = 0; i < count_T - 1; i++)
	{
		sem_destroy(sem_Division + i);
		sem_destroy(sem_Elimination + i);
	}
	sem_destroy(&sem_leader);
	delete[] sem_Division;
	delete[] sem_Elimination;
}

void func(void (*f)())
{
	int counter = 0;
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
		initial(n);
		//平凡算法
		cout << "origin_version: ";
		func(origin_version);
		//动态线程
		cout << "gauss_dynamic: ";
		func(gauss_dynamic);
		//静态线程
		cout << "gauss_static: ";
		func(gauss_static);
		//结合SIMD
		cout << "gauss_SIMD_static: ";
		func(gauss_SIMD_static);
		//列划分
		cout << "gauss_col: ";
		func(gauss_col);
		n *= 2;
		cout << endl;
	}
	return 0;
}
