#include <iostream>
#include <mpi.h>
#include <sys/time.h>
#include <stdlib.h>
#include <pmmintrin.h>
#include <omp.h>
using namespace std;
using namespace MPI;

static const int N = 1000;
static const int task = 1;
static const int thread_count = 4;

float test[N][N] = {0};
float mat[N][N] = {0};

void init_mat(float test[][N])
{
    srand((unsigned)time(NULL));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            test[i][j] = rand() / 100;
}

void reset_mat(float mat[][N], float test[][N])
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            mat[i][j] = test[i][j];
}

void naive_lu(float mat[][N])
{
    for (int k = 0; k < N; k++)
    {
        for (int j = k + 1; j < N; j++)
            mat[k][j] = mat[k][j] / mat[k][k];
        mat[k][k] = 1.0;
        for (int i = k + 1; i < N; i++)
        {
            for (int j = k + 1; j < N; j++)
                mat[i][j] = mat[i][j] - mat[i][k] * mat[k][j];
            mat[i][k] = 0.0;
        }
    }
}

void print_mat(float mat[][N])
{
    if (N > 16)
        return;
    cout << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            cout << mat[i][j] << " ";
        cout << endl;
    }
    cout << endl;
}

class MPI_Block
{
public:
    void eliminate(float mat[][N], int rank, int num_proc)
    {
        int block = N / num_proc;
        //    未能整除划分的剩余部分
        int remain = N % num_proc;

        int begin = rank * block;
        //    当前进程为最后一个进程时，需处理剩余部分
        int end = rank != num_proc - 1 ? begin + block : begin + block + remain;
        for (int k = 0; k < N; k++)
        {
            //        判断当前行是否是自己的任务
            if (k >= begin && k < end)
            {
                for (int j = k + 1; j < N; j++)
                    mat[k][j] = mat[k][j] / mat[k][k];
                mat[k][k] = 1.0;
                //            向之后的进程发送消息
                for (int p = rank + 1; p < num_proc; p++)
                    MPI_Send(mat[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
            }
            else
            {
                int cur_p = k / block;
                //            当所处行属于当前进程前一进程的任务，需接收消息
                if (cur_p < rank)
                    MPI_Recv(mat[k], N, MPI_FLOAT, cur_p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            for (int i = begin; i < end && i < N; i++)
            {
                if (i >= k + 1)
                {
                    for (int j = k + 1; j < N; j++)
                        mat[i][j] = mat[i][j] - mat[i][k] * mat[k][j];
                    mat[i][k] = 0.0;
                }
            }
        }
    }

    void run()
    {
        timeval t_start;
        timeval t_end;

        int num_proc;
        int rank;

        MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int block = N / num_proc;
        int remain = N % num_proc;
        if (rank == 0)
        {
            reset_mat(mat, test);
            gettimeofday(&t_start, NULL);
            //        在0号进程进行任务划分
            for (int i = 1; i < num_proc; i++)
            {
                if (i != num_proc - 1)
                {
                    for (int j = 0; j < block; j++)
                        MPI_Send(mat[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
                }
                else
                {
                    for (int j = 0; j < block + remain; j++)
                        MPI_Send(mat[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
                }
            }
            eliminate(mat, rank, num_proc);
            //        处理完0号进程自己的任务后需接收其他进程处理之后的结果
            for (int i = 1; i < num_proc; i++)
            {
                if (i != num_proc - 1)
                {
                    for (int j = 0; j < block; j++)
                        MPI_Recv(mat[i * block + j], N, MPI_FLOAT, i, 1,
                                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                else
                {
                    for (int j = 0; j < block + remain; j++)
                        MPI_Recv(mat[i * block + j], N, MPI_FLOAT, i, 1,
                                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            gettimeofday(&t_end, NULL);
            cout << "Block MPI LU time cost: "
                 << 1000 * (t_end.tv_sec - t_start.tv_sec) +
                        0.001 * (t_end.tv_usec - t_start.tv_usec)
                 << "ms" << endl;
            print_mat(mat);
        }
        else
        {
            //        非0号进程先接收任务
            if (rank != num_proc - 1)
            {
                for (int j = 0; j < block; j++)
                    MPI_Recv(mat[rank * block + j], N, MPI_FLOAT, 0, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                for (int j = 0; j < block + remain; j++)
                    MPI_Recv(mat[rank * block + j], N, MPI_FLOAT, 0, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            eliminate(mat, rank, num_proc);
            //        处理完后向零号进程返回结果
            if (rank != num_proc - 1)
            {
                for (int j = 0; j < block; j++)
                    MPI_Send(mat[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
            }
            else
            {
                for (int j = 0; j < block + remain; j++)
                    MPI_Send(mat[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
            }
        }
    }

    void eliminate_opt(float mat[][N], int rank, int num_proc)
    {
        __m128 t1, t2, t3;
        int block = N / num_proc;
        int remain = N % num_proc;
        int begin = rank * block;
        int end = rank != num_proc - 1 ? begin + block : begin + block + remain;
#pragma omp parallel num_threads(thread_count)
        for (int k = 0; k < N; k++)
        {
            if (k >= begin && k < end)
            {
                float temp1[4] = {mat[k][k], mat[k][k], mat[k][k], mat[k][k]};
                t1 = _mm_loadu_ps(temp1);
                int j = k + 1;
#pragma omp for schedule(guided, 20)
                for (j; j < N - 3; j += 4)
                {
                    t2 = _mm_loadu_ps(mat[k] + j);
                    t3 = _mm_div_ps(t2, t1);
                    _mm_storeu_ps(mat[k] + j, t3);
                }
#pragma omp for schedule(guided, 20)
                for (j; j < N; j++)
                {
                    mat[k][j] = mat[k][j] / mat[k][k];
                }
                mat[k][k] = 1.0;
                for (int p = rank + 1; p < num_proc; p++)
                    MPI_Send(mat[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
            }
            else
            {
                int cur_p = k / block;
                if (cur_p < rank)
                    MPI_Recv(mat[k], N, MPI_FLOAT, cur_p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            for (int i = begin; i < end && i < N; i++)
            {
                if (i >= k + 1)
                {
                    float temp2[4] = {mat[i][k], mat[i][k], mat[i][k], mat[i][k]};
                    t1 = _mm_loadu_ps(temp2);
                    int j = k + 1;
#pragma omp for schedule(guided, 20)
                    for (j; j <= N - 3; j += 4)
                    {
                        t2 = _mm_loadu_ps(mat[i] + j);
                        t3 = _mm_loadu_ps(mat[k] + j);
                        t3 = _mm_mul_ps(t1, t3);
                        t2 = _mm_sub_ps(t2, t3);
                        _mm_storeu_ps(mat[i] + j, t2);
                    }
#pragma omp for schedule(guided, 20)
                    for (j; j < N; j++)
                        mat[i][j] = mat[i][j] - mat[i][k] * mat[k][j];
                    mat[i][k] = 0;
                }
            }
        }
    }

    void run_opt()
    {
        timeval t_start;
        timeval t_end;

        int num_proc;
        int rank;

        MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int block = N / num_proc;
        int remain = N % num_proc;
        if (rank == 0)
        {
            reset_mat(mat, test);
            gettimeofday(&t_start, NULL);
            for (int i = 1; i < num_proc; i++)
            {
                if (i != num_proc - 1)
                {
                    for (int j = 0; j < block; j++)
                        MPI_Send(mat[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
                }
                else
                {
                    for (int j = 0; j < block + remain; j++)
                        MPI_Send(mat[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
                }
            }
            eliminate_opt(mat, rank, num_proc);
            for (int i = 1; i < num_proc; i++)
            {
                if (i != num_proc - 1)
                {
                    for (int j = 0; j < block; j++)
                        MPI_Recv(mat[i * block + j], N, MPI_FLOAT, i, 1,
                                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                else
                {
                    for (int j = 0; j < block + remain; j++)
                        MPI_Recv(mat[i * block + j], N, MPI_FLOAT, i, 1,
                                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            gettimeofday(&t_end, NULL);
            cout << "Block MPI LU with SSE and OpenMP time cost: "
                 << 1000 * (t_end.tv_sec - t_start.tv_sec) +
                        0.001 * (t_end.tv_usec - t_start.tv_usec)
                 << "ms" << endl;
            print_mat(mat);
        }
        else
        {
            if (rank != num_proc - 1)
            {
                for (int j = 0; j < block; j++)
                    MPI_Recv(mat[rank * block + j], N, MPI_FLOAT, 0, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                for (int j = 0; j < block + remain; j++)
                    MPI_Recv(mat[rank * block + j], N, MPI_FLOAT, 0, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            eliminate_opt(mat, rank, num_proc);
            if (rank != num_proc - 1)
            {
                for (int j = 0; j < block; j++)
                    MPI_Send(mat[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
            }
            else
            {
                for (int j = 0; j < block + remain; j++)
                    MPI_Send(mat[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
            }
        }
    }
};

class MPI_Recycle
{
public:
    void eliminate(float mat[][N], int rank, int num_proc)
    {
        //    所有进程进行1次迭代的计算行数
        int seg = task * num_proc;
        for (int k = 0; k < N; k++)
        {
            //        判断当前行是否是自己的任务
            if (int((k % seg) / task) == rank)
            {
                for (int j = k + 1; j < N; j++)
                    mat[k][j] = mat[k][j] / mat[k][k];
                mat[k][k] = 1.0;
                //            完成计算后向其他进程发送消息
                for (int p = 0; p < num_proc; p++)
                    if (p != rank)
                        MPI_Send(mat[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
            }
            else
            {
                //            如果当前行不是自己的任务，接收来自当前行处理进程的消息
                MPI_Recv(mat[k], N, MPI_FLOAT, int((k % seg) / task), 2,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            for (int i = k + 1; i < N; i++)
            {
                if (int((i % seg) / task) == rank)
                {
                    for (int j = k + 1; j < N; j++)
                        mat[i][j] = mat[i][j] - mat[i][k] * mat[k][j];
                    mat[i][k] = 0.0;
                }
            }
        }
    }

    void run()
    {
        timeval t_start;
        timeval t_end;

        int num_proc;
        int rank;

        MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int seg = task * num_proc;
        if (rank == 0)
        {
            reset_mat(mat, test);
            gettimeofday(&t_start, NULL);
            //        在0号进程进行任务划分
            for (int i = 0; i < N; i++)
            {
                int flag = (i % seg) / task;
                if (flag == rank)
                    continue;
                else
                    MPI_Send(mat[i], N, MPI_FLOAT, flag, 0, MPI_COMM_WORLD);
            }
            eliminate(mat, rank, num_proc);
            //        处理完0号进程自己的任务后需接收其他进程处理之后的结果
            for (int i = 0; i < N; i++)
            {
                int flag = (i % seg) / task;
                if (flag == rank)
                    continue;
                else
                    MPI_Recv(mat[i], N, MPI_FLOAT, flag, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            gettimeofday(&t_end, NULL);
            cout << "Recycle MPI LU time cost: "
                 << 1000 * (t_end.tv_sec - t_start.tv_sec) +
                        0.001 * (t_end.tv_usec - t_start.tv_usec)
                 << "ms" << endl;
            print_mat(mat);
        }
        else
        {
            //        非0号进程先接收任务
            for (int i = task * rank; i < N; i += seg)
            {
                for (int j = 0; j < task && i + j < N; j++)
                    MPI_Recv(mat[i + j], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            eliminate(mat, rank, num_proc);
            //        处理完后向零号进程返回结果
            for (int i = task * rank; i < N; i += seg)
            {
                for (int j = 0; j < task && i + j < N; j++)
                    MPI_Send(mat[i + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
            }
        }
    }

    void eliminate_opt(float mat[][N], int rank, int num_proc)
    {
        __m128 t1, t2, t3;
        int seg = task * num_proc;
#pragma omp parallel num_threads(thread_count)
        for (int k = 0; k < N; k++)
        {
            if (int((k % seg) / task) == rank)
            {
                float temp1[4] = {mat[k][k], mat[k][k], mat[k][k], mat[k][k]};
                t1 = _mm_loadu_ps(temp1);
                int j = k + 1;
#pragma omp for schedule(guided, 20)
                for (j; j < N - 3; j += 4)
                {
                    t2 = _mm_loadu_ps(mat[k] + j);
                    t3 = _mm_div_ps(t2, t1);
                    _mm_storeu_ps(mat[k] + j, t3);
                }
#pragma omp for schedule(guided, 20)
                for (j; j < N; j++)
                {
                    mat[k][j] = mat[k][j] / mat[k][k];
                }
                mat[k][k] = 1.0;

                for (int p = 0; p < num_proc; p++)
                    if (p != rank)
                        MPI_Send(mat[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Recv(mat[k], N, MPI_FLOAT, int((k % seg) / task), 2,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            for (int i = k + 1; i < N; i++)
            {
                if (int((i % seg) / task) == rank)
                {
                    float temp2[4] = {mat[i][k], mat[i][k], mat[i][k], mat[i][k]};
                    t1 = _mm_loadu_ps(temp2);
                    int j = k + 1;
#pragma omp for schedule(guided, 20)
                    for (j; j <= N - 3; j += 4)
                    {
                        t2 = _mm_loadu_ps(mat[i] + j);
                        t3 = _mm_loadu_ps(mat[k] + j);
                        t3 = _mm_mul_ps(t1, t3);
                        t2 = _mm_sub_ps(t2, t3);
                        _mm_storeu_ps(mat[i] + j, t2);
                    }
#pragma omp for schedule(guided, 20)
                    for (j; j < N; j++)
                        mat[i][j] = mat[i][j] - mat[i][k] * mat[k][j];
                    mat[i][k] = 0;
                }
            }
        }
    }

    void run_opt()
    {
        timeval t_start;
        timeval t_end;

        int num_proc;
        int rank;

        MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int seg = task * num_proc;
        if (rank == 0)
        {
            reset_mat(mat, test);
            gettimeofday(&t_start, NULL);
            for (int i = 0; i < N; i++)
            {
                int flag = (i % seg) / task;
                if (flag == rank)
                    continue;
                else
                    MPI_Send(mat[i], N, MPI_FLOAT, flag, 0, MPI_COMM_WORLD);
            }
            eliminate_opt(mat, rank, num_proc);
            for (int i = 0; i < N; i++)
            {
                int flag = (i % seg) / task;
                if (flag == rank)
                    continue;
                else
                    MPI_Recv(mat[i], N, MPI_FLOAT, flag, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            gettimeofday(&t_end, NULL);
            cout << "Recycle MPI LU with SSE and OpenMP time cost: "
                 << 1000 * (t_end.tv_sec - t_start.tv_sec) +
                        0.001 * (t_end.tv_usec - t_start.tv_usec)
                 << "ms" << endl;
            print_mat(mat);
        }
        else
        {
            for (int i = task * rank; i < N; i += seg)
            {
                for (int j = 0; j < task && i + j < N; j++)
                    MPI_Recv(mat[i + j], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            eliminate_opt(mat, rank, num_proc);
            for (int i = task * rank; i < N; i += seg)
            {
                for (int j = 0; j < task && i + j < N; j++)
                    MPI_Send(mat[i + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
            }
        }
    }
};

class MPI_Pipeline
{
public:
    void eliminate(float mat[][N], int rank, int num_proc)
    {
        int seg = task * num_proc;
        //    计算当前进程的前一进程及下一进程
        int pre_proc = (rank + (num_proc - 1)) % num_proc;
        int next_proc = (rank + 1) % num_proc;
        for (int k = 0; k < N; k++)
        {
            //        判断当前行是否是自己的任务
            if (int((k % seg) / task) == rank)
            {
                for (int j = k + 1; j < N; j++)
                    mat[k][j] = mat[k][j] / mat[k][k];
                mat[k][k] = 1.0;
                //            处理完自己的任务后向下一进程发送消息
                MPI_Send(mat[k], N, MPI_FLOAT, next_proc, 2, MPI_COMM_WORLD);
            }
            else
            {
                //            如果当前行不是当前进程的任务，则接收前一进程的消息
                MPI_Recv(mat[k], N, MPI_FLOAT, pre_proc, 2,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //            如果当前行不是下一进程的任务，需将消息进行传递
                if (int((k % seg) / task) != next_proc)
                    MPI_Send(mat[k], N, MPI_FLOAT, next_proc, 2, MPI_COMM_WORLD);
            }
            for (int i = k + 1; i < N; i++)
            {
                if (int((i % seg) / task) == rank)
                {
                    for (int j = k + 1; j < N; j++)
                        mat[i][j] = mat[i][j] - mat[i][k] * mat[k][j];
                    mat[i][k] = 0.0;
                }
            }
        }
    }

    void run()
    {
        timeval t_start;
        timeval t_end;

        int num_proc;
        int rank;

        MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int seg = task * num_proc;
        if (rank == 0)
        {
            reset_mat(mat, test);
            gettimeofday(&t_start, NULL);
            //        在0号进程进行任务划分
            for (int i = 0; i < N; i++)
            {
                int flag = (i % seg) / task;
                if (flag == rank)
                    continue;
                else
                    MPI_Send(mat[i], N, MPI_FLOAT, flag, 0, MPI_COMM_WORLD);
            }
            eliminate(mat, rank, num_proc);
            //        处理完0号进程自己的任务后需接收其他进程处理之后的结果
            for (int i = 0; i < N; i++)
            {
                int flag = (i % seg) / task;
                if (flag == rank)
                    continue;
                else
                    MPI_Recv(mat[i], N, MPI_FLOAT, flag, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            gettimeofday(&t_end, NULL);
            cout << "Pipeline MPI LU time cost: "
                 << 1000 * (t_end.tv_sec - t_start.tv_sec) +
                        0.001 * (t_end.tv_usec - t_start.tv_usec)
                 << "ms" << endl;
            print_mat(mat);
        }
        else
        {
            //        非0号进程先接收任务
            for (int i = task * rank; i < N; i += seg)
            {
                for (int j = 0; j < task && i + j < N; j++)
                    MPI_Recv(mat[i + j], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            eliminate(mat, rank, num_proc);
            //        处理完后向零号进程返回结果
            for (int i = task * rank; i < N; i += seg)
            {
                for (int j = 0; j < task && i + j < N; j++)
                    MPI_Send(mat[i + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
            }
        }
    }

    void eliminate_opt(float mat[][N], int rank, int num_proc)
    {
        __m128 t1, t2, t3;
        int seg = task * num_proc;
        int pre_proc = (rank + (num_proc - 1)) % num_proc;
        int next_proc = (rank + 1) % num_proc;
#pragma omp parallel num_threads(thread_count)
        for (int k = 0; k < N; k++)
        {
            if (int((k % seg) / task) == rank)
            {
                float temp1[4] = {mat[k][k], mat[k][k], mat[k][k], mat[k][k]};
                t1 = _mm_loadu_ps(temp1);
                int j = k + 1;
#pragma omp for schedule(guided, 20)
                for (j; j < N - 3; j += 4)
                {
                    t2 = _mm_loadu_ps(mat[k] + j);
                    t3 = _mm_div_ps(t2, t1);
                    _mm_storeu_ps(mat[k] + j, t3);
                }
#pragma omp for schedule(guided, 20)
                for (j; j < N; j++)
                {
                    mat[k][j] = mat[k][j] / mat[k][k];
                }
                mat[k][k] = 1.0;
                MPI_Send(mat[k], N, MPI_FLOAT, next_proc, 2, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Recv(mat[k], N, MPI_FLOAT, pre_proc, 2,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (int((k % seg) / task) != next_proc)
                    MPI_Send(mat[k], N, MPI_FLOAT, next_proc, 2, MPI_COMM_WORLD);
            }
            for (int i = k + 1; i < N; i++)
            {
                if (int((i % seg) / task) == rank)
                {
                    float temp2[4] = {mat[i][k], mat[i][k], mat[i][k], mat[i][k]};
                    t1 = _mm_loadu_ps(temp2);
                    int j = k + 1;
#pragma omp for schedule(guided, 20)
                    for (j; j <= N - 3; j += 4)
                    {
                        t2 = _mm_loadu_ps(mat[i] + j);
                        t3 = _mm_loadu_ps(mat[k] + j);
                        t3 = _mm_mul_ps(t1, t3);
                        t2 = _mm_sub_ps(t2, t3);
                        _mm_storeu_ps(mat[i] + j, t2);
                    }
#pragma omp for schedule(guided, 20)
                    for (j; j < N; j++)
                        mat[i][j] = mat[i][j] - mat[i][k] * mat[k][j];
                    mat[i][k] = 0;
                }
            }
        }
    }

    void run_opt()
    {
        timeval t_start;
        timeval t_end;

        int num_proc;
        int rank;

        MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int seg = task * num_proc;
        if (rank == 0)
        {
            reset_mat(mat, test);
            gettimeofday(&t_start, NULL);
            for (int i = 0; i < N; i++)
            {
                int flag = (i % seg) / task;
                if (flag == rank)
                    continue;
                else
                    MPI_Send(mat[i], N, MPI_FLOAT, flag, 0, MPI_COMM_WORLD);
            }
            eliminate_opt(mat, rank, num_proc);
            for (int i = 0; i < N; i++)
            {
                int flag = (i % seg) / task;
                if (flag == rank)
                    continue;
                else
                    MPI_Recv(mat[i], N, MPI_FLOAT, flag, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            gettimeofday(&t_end, NULL);
            cout << "Pipeline MPI LU with SSE and OpenMP time cost: "
                 << 1000 * (t_end.tv_sec - t_start.tv_sec) +
                        0.001 * (t_end.tv_usec - t_start.tv_usec)
                 << "ms" << endl;
            print_mat(mat);
        }
        else
        {
            for (int i = task * rank; i < N; i += seg)
            {
                for (int j = 0; j < task && i + j < N; j++)
                    MPI_Recv(mat[i + j], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            eliminate_opt(mat, rank, num_proc);
            for (int i = task * rank; i < N; i += seg)
            {
                for (int j = 0; j < task && i + j < N; j++)
                    MPI_Send(mat[i + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
            }
        }
    }
};

int main()
{
    init_mat(test);

    MPI_Init(NULL, NULL);

    MPI_Block mpi_block;
    mpi_block.run();
    mpi_block.run_opt();

    MPI_Recycle mpi_recycle;
    mpi_recycle.run();
    mpi_recycle.run_opt();

    MPI_Pipeline mpi_pipeline;
    mpi_pipeline.run();
    mpi_pipeline.run_opt();

    MPI_Finalize();
    return 0;
}