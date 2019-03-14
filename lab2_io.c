#include "lab2_io.h"

void read_matrix (const char* input_filename, int* M, int* N, float** D){
	FILE *fin = fopen(input_filename, "r");

	fscanf(fin, "%d%d", M, N);
	
	int num_elements = (*M) * (*N);
	*D = (float*) malloc(sizeof(float)*(num_elements));
	
	for (int i = 0; i < num_elements; i++){
		fscanf(fin, "%f", (*D + i));
	}
	fclose(fin);
}

void write_result (int M, 
		int N, 
		float* D, 
		float* U, 
		float* SIGMA, 
		float* V_T,
		int K, 
		float* D_HAT,
		double computation_time){
	// Will contain output code
}
