#include "pch.h"
#include <iostream>
#include "mpi.h"
#include <iomanip>
#include <ctime>


using namespace std;

int n;
int row;
int root = 0, rem_row, tmp;
double** matrix;
double* b;
double* solution;
int chunk_row;

void Create_Matrix(int m_size)
{
	srand(time(NULL));
	n = m_size;
	matrix = new double*[n];

	for (int i = 0; i < n; i++)
	{
		matrix[i] = new double[n];
		memset(matrix[i], false, n * sizeof(double));
	}


	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			matrix[i][j] = 1 + rand() % 1000;
		}
	}
}

void Create_b_Vector(int m_size)
{
	srand(time(NULL));
	n = m_size;
	b = (double*)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++)
	{
		b[i] = 1 + rand() % 1000;
	}
}

void Print_Matrix()
{

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf("%.2f", matrix[i][j]);
			printf("\t");
		}
		printf("%.2f", b[i]);
		printf("\n");
	}
	printf("\n");
}

int main(int argc, char* argv[])
{
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int chunk;
	double pivot_b;
	int pivot = 0;

	n = 3 + rand() % 7;

	int* b_chnk_vec;
	b_chnk_vec = (int*)malloc(size * sizeof(int));

	int* matrix_disp_vec;
	matrix_disp_vec = (int*)malloc(size * sizeof(int));


	int* matrix_chnk_vec;
	matrix_chnk_vec = (int*)malloc(size * sizeof(int));

	int* b_disp_vec;
	b_disp_vec = (int*)malloc(size * sizeof(int));

	double* local_vector;

	clock_t duration = clock();

	while ((n - pivot - 1) / size > 0)
	{
		double* pivot_eqn;
		pivot_eqn = (double*)malloc((n - pivot) * sizeof(double));

		double* sendbuf;
		sendbuf = (double*)malloc((n - pivot - 1) * (n - pivot) * sizeof(double));

		if (rank == 0 && pivot == 0)
		{
			Create_Matrix(n);
			Create_b_Vector(n);
			Print_Matrix();
		}

		if (rank == 0)	// master operations
		{
			int k = 0;

			for (int i = 0; i < n - pivot - 1; i++)
			{
				for (int j = 0; j < n - pivot; j++)
				{
					sendbuf[k] = matrix[pivot + i + 1][pivot + j];
					k++;
				}
			}

			if ((n - pivot - 1) % size != 0)
			{
				row = round((n - pivot - 1) / size);
				rem_row = (n - pivot - 1) - row * (size - 1);
				tmp = 0;
				for (int i = 0; i < size; i++)
				{
					if (i != size - 1)
					{
						matrix_chnk_vec[i] = row * (n - pivot);
						b_chnk_vec[i] = row;
					}
					else
					{
						matrix_chnk_vec[i] = rem_row * (n - pivot);
						b_chnk_vec[i] = rem_row;
						rem_row = 0;
					}
				}
			}
			else
			{
				row = (n - pivot - 1) / size;
				for (int i = 0; i < size; i++)
				{
					matrix_chnk_vec[i] = row * (n - pivot);
					b_chnk_vec[i] = row;
				}
			}

			// calculate dispalcement vectors

			matrix_disp_vec[0] = 0;
			b_disp_vec[0] = pivot + 1;
			for (int i = 1; i < size; i++)
			{
				matrix_disp_vec[i] = matrix_disp_vec[i - 1] + matrix_chnk_vec[i - 1];
				b_disp_vec[i] = b_disp_vec[i - 1] + b_chnk_vec[i - 1];
			}

			//fetch the pivot equation
			for (int j = 0; j < n - pivot; j++)
			{
				pivot_eqn[j] = matrix[pivot][j + pivot];
			}

			pivot_b = b[pivot];
		}

		MPI_Bcast(&n, 1, MPI_INT, root, MPI_COMM_WORLD);

		//distribute the pivot equation

		MPI_Bcast(pivot_eqn, n - pivot, MPI_DOUBLE, root, MPI_COMM_WORLD);
		MPI_Bcast(&pivot_b, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

		//distribute the chunk vector
		MPI_Scatter(matrix_chnk_vec, 1, MPI_INTEGER, &chunk, 1, MPI_INT, root, MPI_COMM_WORLD);

		chunk_row = n - pivot;

		double* local_b;
		local_b = (double*)malloc(chunk / chunk_row * sizeof(double));

		double** local_matrix;


		local_vector = (double*)malloc(chunk * sizeof(double));

		local_matrix = new double *[chunk / chunk_row];
		for (int i = 0; i < chunk / chunk_row; i++)
		{
			local_matrix[i] = new double[chunk_row];
			memset(local_matrix[i], false, chunk_row * sizeof(double));
		}


		//distribute the data
		MPI_Scatterv(sendbuf, matrix_chnk_vec, matrix_disp_vec, MPI_DOUBLE, local_vector, chunk, MPI_DOUBLE, root, MPI_COMM_WORLD);
		MPI_Scatterv(b, b_chnk_vec, b_disp_vec, MPI_DOUBLE, local_b, chunk / chunk_row, MPI_DOUBLE, root, MPI_COMM_WORLD);

		int k = 0;

		for (int i = 0; i < chunk / chunk_row; i++)
		{
			for (int j = 0; j < chunk_row; j++)
			{
				local_matrix[i][j] = local_vector[k];
				k++;
			}
		}


		double t;

		for (int i = 0; i < chunk / chunk_row; i++)
		{
			t = local_matrix[i][0] / pivot_eqn[0];
			for (int j = 0; j < chunk_row; j++)
			{
				local_matrix[i][j] -= (t * pivot_eqn[j]);

			}
			local_b[i] -= (t * pivot_b);
		}


		k = 0;

		for (int i = 0; i < chunk / chunk_row; i++)
		{
			for (int j = 0; j < chunk_row; j++)
			{
				local_vector[k] = local_matrix[i][j];
				k++;
			}
		}

		MPI_Gatherv(local_vector, chunk, MPI_DOUBLE, sendbuf, matrix_chnk_vec, matrix_disp_vec, MPI_DOUBLE, root, MPI_COMM_WORLD);
		MPI_Gatherv(local_b, chunk / chunk_row, MPI_DOUBLE, b, b_chnk_vec, b_disp_vec, MPI_DOUBLE, root, MPI_COMM_WORLD);

		if (rank == 0)
		{
			int k = 0;

			for (int i = 0; i < n - pivot - 1; i++)
			{
				for (int j = 0; j < n - pivot; j++)
				{
					matrix[pivot + i + 1][pivot + j] = sendbuf[k];
					k++;
				}
			}

		}

		pivot++;

		MPI_Barrier(MPI_COMM_WORLD);

	}

	MPI_Finalize();

	double t;

	for (int k = pivot; k < n - 1; k++)
	{
		for (int i = k + 1; i < n; i++)
		{
			t = matrix[i][k] / matrix[k][k];
			for (int j = k; j < n; j++)
			{
				matrix[i][j] -= (t * matrix[k][j]);
			}
			b[i] -= (t * b[k]);
		}
	}

	Print_Matrix();

	duration = clock() - duration;

	solution = (double*)malloc(n * sizeof(double));

	solution[n - 1] = b[n - 1] / matrix[n - 1][n - 1];

	for (int i = n - 2; i >= 0; i--)
	{
		solution[i] = b[n - 1];
		for (int j = i + 1; j < n; j++)

			solution[i] = solution[i] - matrix[i][j] * solution[j];
		solution[i] = solution[i] / matrix[i][i];
	}
	for (int i = 0; i < n; i++)
	{
		printf("x[%d] = %.2f", i, solution[i]);
		printf("\n");
	}

	printf("Duration of the Proccess is %d mseconds", duration);


	return 0;
}
