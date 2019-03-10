#include "lab2_io.h"
#include "lab2_omp.h"

#include <stdlib.h>
#include <omp.h>

/*
	Arguments:
		arg1: input filename (consist M, N and D)
		arg2: retention (percentage of information to be retained by PCA)
*/

int main(int argc, char const *argv[])
{
	if (argc < 3){
		printf("\nLess Arguments\n");
		return 0;
	}

	if (argc > 3){
		printf("\nTOO many Arguments\n");
		return 0;
	}

	//---------------------------------------------------------------------
	int M;			//no of rows (samples) in input matrix D (input)
	int N;			//no of columns (features) in input matrix D (input)
	float* D;		//1D array of M x N matrix to be reduced (input)
	float* U;		//1D array of N x N matrix U (to be computed by SVD)
	float* SIGMA;	//1D array of N x M diagonal matrix SIGMA (to be computed by SVD)
	float* V_T;		//1D array of M x M matrix V_T (to be computed by SVD)
	int K;			//no of coulmns (features) in reduced matrix D_HAT (to be computed by PCA)
	float *D_HAT;	//1D array of M x K reduced matrix (to be computed by PCA)
	int retention;	//percentage of information to be retained by PCA (command line input)
	//---------------------------------------------------------------------

	retention = atoi(argv[2]);	//retention = 90 means 90% of information should be retained

	float start_time, end_time;
	double computation_time;

	/*
		-- Pre-defined function --
		reads matrix and its dimentions from input file and creats array D
	    #elements in D is M * N
        format - 
        --------------------------------------------------------------------------------------
        | D[0][0] | D[0][1] | ... | D[0][N-1] | D[1][0] | ... | D[1][N-1] | ... | D[M-1][N-1] |
        --------------------------------------------------------------------------------------
	*/
	read_matrix (argv[1], &M, &N, &D);

	U = (float*) malloc(sizeof(float) * N*N);
	SIGMA = (float*) malloc(sizeof(float) * N);
	V_T = (float*) malloc(sizeof(float) * M*M);

	start_time = omp_get_wtime();
	
	// /*
	// 	*****************************************************
	// 		TODO -- You must implement these two function
	// 	*****************************************************
	// */
	SVD(M, N, D, &U, &SIGMA, &V_T);
	PCA(retention, M, N, D, U, SIGMA, &D_HAT, &K);

	end_time = omp_get_wtime();
	computation_time = ((double) (end_time - start_time));
	
	/*
		--Pre-defined functions --
		checks for correctness of results computed by SVD and PCA
		and outputs the results
	*/
	write_result(M, N, D, U, SIGMA, V_T, K, D_HAT, computation_time);

	return 0;
}
