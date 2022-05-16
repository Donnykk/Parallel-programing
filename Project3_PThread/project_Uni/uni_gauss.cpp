#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>
#include <arm_neon.h>
#include <sys/time.h>
#define millitime(x) (x.tv_sec * 1000 + x.tv_usec / 1000.0)
#define num_col 130
#define num_row 22
#define num_ele 8
using namespace std;

const int totalCols = num_col / 32 + (num_col % 32 != 0 ? 1 : 0);
const int totalRows = num_row + num_ele;
int **bitmap = NULL;
int *lp = NULL;
int count_T = 4;
map<int, int> lp_Row;

typedef struct
{
	int t_id; //线程id
} threadParam_t;

void readFile(string path, int &currentRow)
{
	string data;
	stringstream ss;
	ss.clear();
	ifstream input(path);
	while (getline(input, data))
	{
		ss << data;
		int index;
		while (ss >> index)
		{
			bitmap[currentRow][index / 32] |= (1 << (index % 32));
			lp[currentRow] = max(index, lp[currentRow]);
		}
		if (currentRow >= num_ele)
			lp_Row[lp[currentRow]] = currentRow;
		currentRow++;
		ss.clear();
	}
	input.close();
}

void init()
{
	if (bitmap != NULL)
	{
		for (int i = 0; i < totalRows; i++)
			delete[] bitmap[i];
		delete[] bitmap;
	}
	if (lp != NULL)
		delete[] lp;
	bitmap = new int *[totalRows];
	for (int i = 0; i < totalRows; i++)
		bitmap[i] = new int[totalCols]{0};
	lp = new int[totalRows]{-1};
	lp_Row.clear();
	int currentRow = 0;
	readFile("2.txt", currentRow);
	readFile("1.txt", currentRow);
}

void Uni_Gauss()
{
	for (int i = 0; i < num_ele; i++)
	{
		while (lp[i] > -1)
		{
			if (!(lp_Row.find(lp[i]) == lp_Row.end()))
			{
				int rowR = lp_Row[lp[i]];
				bool modified = false;
				for (int j = 0; j < totalCols - 1; j++)
					bitmap[i][j] ^= bitmap[rowR][j];
				for (int j = totalCols - 1; j >= 0; j--)
				{
					if (bitmap[i][j] != 0)
					{
						modified = true;
						for (int k = 31; k >= 0; k--)
						{
							if ((bitmap[i][j] & (1 << k)) != 0)
							{
								lp[i] = j * 32 + k;
								break;
							}
						}
						break;
					}
				}
				if (modified == false)
					lp[i] = -1;
			}
			else
			{
				lp_Row[lp[i]] = i;
				break;
			}
		}
	}
}

pthread_barrier_t barrier1;
pthread_barrier_t barrier2;
pthread_barrier_t barrier3;

void *threadFunc(void *param)
{
	long long t_id = (long long)param;
	for (int i = 0; i < num_ele; i++)
	{
		while (lp[i] > -1)
		{
			if (!(lp_Row.find(lp[i]) == lp_Row.end()))
			{
				int rowR = lp_Row[lp[i]];
				for (int j = t_id; j < totalCols - 1; j += count_T)
					bitmap[i][j] ^= bitmap[rowR][j];
			}
			else
			{
				pthread_barrier_wait(&barrier3);
				lp_Row[lp[i]] = i;
				break;
			}
			pthread_barrier_wait(&barrier1);
			if (t_id == 0)
			{
				bool modified = false;
				for (int j = totalCols - 1; j >= 0; j--)
				{
					if (bitmap[i][j] != 0)
					{
						modified = true;
						for (int k = 31; k >= 0; k--)
						{
							if ((bitmap[i][j] & (1 << k)) != 0)
							{
								lp[i] = j * 32 + k;
								break;
							}
						}
						break;
					}
				}
				if (modified == false)
					lp[i] = -1;
			}
			pthread_barrier_wait(&barrier2);
		}
	}
	pthread_exit(NULL);
}

void PThread_Uni_Gauss()
{
	pthread_barrier_init(&barrier1, NULL, count_T);
	pthread_barrier_init(&barrier2, NULL, count_T);
	pthread_barrier_init(&barrier3, NULL, count_T);
	pthread_t handles[count_T];
	threadParam_t param[count_T];
	for (int t_id = 0; t_id < count_T; t_id++)
	{
		param[t_id].t_id = t_id;
		pthread_create(handles + t_id, NULL, threadFunc, (void *)(long long)t_id);
	}
	for (int i = 0; i < count_T; i++)
		pthread_join(handles[i], NULL);
	pthread_barrier_destroy(&barrier1);
	pthread_barrier_destroy(&barrier2);
	pthread_barrier_destroy(&barrier3);
}

void *threadFunc_plus(void *param)
{
	long long t_id = (long long)param;
	for (int i = 0; i < num_ele; i++)
	{
		while (lp[i] > -1)
		{
			int j = t_id;
			if (!(lp_Row.find(lp[i]) == lp_Row.end()))
			{
				int rowR = lp_Row[lp[i]];
				int32x4_t v0, v1;
				for (; j <= totalCols - 4 * count_T; j += 4 * count_T)
				{
					v0 = vld1q_s32(bitmap[i] + j);
					v1 = vld1q_s32(bitmap[rowR] + j);
					v0 = veorq_s32(v0, v1);
					vst1q_s32(bitmap[i] + j, v0);
				}
				for (; j < totalCols - 1; j += count_T)
					bitmap[i][j] ^= bitmap[rowR][j];
			}
			else
			{
				pthread_barrier_wait(&barrier3);
				lp_Row[lp[i]] = i;
				break;
			}
			pthread_barrier_wait(&barrier1);
			if (t_id == 0)
			{
				bool modified = false;
				for (j = totalCols - 1; j >= 3; j -= 4)
				{
					int x = ((bitmap[i][j] | bitmap[i][j - 1]) | (bitmap[i][j - 2] | bitmap[i][j - 3]));
					if (x == 0)
						continue;
					break;
				}
				for (; j >= 0; j--)
				{
					if (bitmap[i][j] != 0)
					{
						modified = true;
						for (int k = 31; k >= 0; k--)
						{
							if ((bitmap[i][j] & (1 << k)) != 0)
							{
								lp[i] = j * 32 + k;
								break;
							}
						}
						break;
					}
				}
				if (modified == false)
					lp[i] = -1;
			}
			pthread_barrier_wait(&barrier2);
		}
	}
	pthread_exit(NULL);
}

void PThread_Neon_Uni_Gauss()
{
	pthread_barrier_init(&barrier1, NULL, count_T);
	pthread_barrier_init(&barrier2, NULL, count_T);
	pthread_barrier_init(&barrier3, NULL, count_T);
	pthread_t handles[count_T];
	threadParam_t param[count_T];
	for (int t_id = 0; t_id < count_T; t_id++)
	{
		param[t_id].t_id = t_id;
		pthread_create(handles + t_id, NULL, threadFunc_plus, (void *)(long long)t_id);
	}
	for (int i = 0; i < count_T; i++)
		pthread_join(handles[i], NULL);
	pthread_barrier_destroy(&barrier1);
	pthread_barrier_destroy(&barrier2);
	pthread_barrier_destroy(&barrier3);
}

void check()
{
	ifstream input("3.txt");
	string data;
	stringstream ss;
	for (int i = 0; i < num_ele; i++)
	{
		getline(input, data);
		if (lp[i] == -1)
		{
			if (data.length() > 1)
			{
				cout << "WRONG";
				return;
			}
			continue;
		}
		ss.clear();
		ss << data;
		int index;
		for (int j = totalCols - 1; j >= 0; j--)
		{
			for (int k = 31; k >= 0; k--)
			{
				if ((bitmap[i][j] & (1 << k)) != 0)
				{
					ss >> index;
					if (index != j * 32 + k)
					{
						cout << "WRONG" << endl;
						return;
					}
				}
			}
		}
	}
	cout << "SUCCESS" << endl;
	input.close();
}

void timer(void (*func)())
{
	int counter = 0;
	timeval begin, start, finish, now;
	float milliseconds = 0, single_time;
	gettimeofday(&begin, NULL);
	gettimeofday(&now, NULL);
	while (millitime(now) - millitime(begin) < 20)
	{
		counter++;
		init();
		gettimeofday(&start, NULL);
		func();
		gettimeofday(&finish, NULL);
		milliseconds += (millitime(finish) - millitime(start));
		gettimeofday(&now, NULL);
	}
	single_time = milliseconds / counter;
	cout << "运行次数：" << counter << " 平均时间：" << single_time << endl;
}

int main()
{
	cout << "线程数:4" <<endl;
	timer(Uni_Gauss);
	check();
	timer(PThread_Uni_Gauss);
	check();
	timer(PThread_Neon_Uni_Gauss);
	check();
	count_T = 8;
	cout << "线程数:8" <<endl;
	timer(Uni_Gauss);
	check();
	timer(PThread_Uni_Gauss);
	check();
	timer(PThread_Neon_Uni_Gauss);
	check();
	return 0;
}
