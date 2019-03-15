#include <malloc.h>
#include <omp.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define TOLERANCE 1e-3
// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */

void multiply(float *M_1, int m1, int n1, float *M_2, int m2, int n2, float *result)
{
    assert(n1 == m2);
    float sum = 0.0;
    int i,j,k;
    for (i=0; i<m1; i++)
    {
        for (j=0; j<n2; j++)
        {
            for (k=0; k<n1; k++)
            {
                sum+=M_1[i*n1+k]*M_2[k*n2+j];
            }
            result[i*n2+j] = sum;
            sum=0.0;
        }        
    }
}

void transpose(float* M, int m, int n, float* M_T)
{
    int i,j;
    for (i=0; i<m; i++)
    {
        for (j=0; j<n; j++)
        {
            M_T[j*m+i]=M[i*m+j];
        }
    }
}

float* initialize_identity(int size)
{
    float *E_0 = (float *) calloc(P*P, sizeof(float));
    for (int i=0; i<size; i++)
    {
        E_0[i*size+i] = 1.0;
    }
}

float l2_norm(float *v_col, int length)
{
    float norm = 0.0, sum_sq=0.0;
    for (int i=0; i<length; i++)
    {
        sum_sq+=v_col[i]*v_col[i];
    }
    return norm = sqrtf(sum_sq);
}

float l2_norm_diagonal_diff(float *A_next, float *A_current, int P)
{
    float norm, sum_sq=0.0;
    for (int i=0; i<P; i++)
    {
        sum_sq+=(A_next[i*P+i]-A_current[i*P+i])*(A_next[i*P+i]-A_current[i*P+i]);
    }
    return norm=sqrtf(sum_sq);
}

void classicalGS(float *A_current, float* A_T, int P, float *Q_current, float *R_current)
{
    float *v_col = (float *) malloc(sizeof(float)*P);
    int col, row, row_;
    float result;
    for (col=0; col<P; col++)
    {
        memcpy(v_col, A_T+col*P, P*sizeof(float));
        //v_col = a_col
        for (row=0; row<col; row++)
        {
            //perform inner product:
            result = 0.0;
            for (row_=0; row_<P; row_++)
            {
                result+=(Q_current[row_*P+row]*(A_T[col*P+row_]));
            }
            R_current[row*P+col]=result;
        }
        R_current[col*P+col] = l2_norm(v_col, P);
        for (row=0; row<P; row++)
        {
            Q_current[row*P+col]=v_col[row]/R_current[col*P+col];
        }
    }
    free(v_col);
}

void SVD(int N, int P, float* D, float** U, float** SIGMA, float** V_T)
{
    /* 1.Perform SVD for D_T */
    // Get eigen-values & eigen-vectors for D_T*D
    float *D_T = (float *) malloc(sizeof(float)*N*P);   
    transpose(D, N, P, D_T);
    float *A = (float *) calloc(P*P, sizeof(float)); //A=D_T*D|(PxP)
    float *A_T = (float *) calloc(P*P, sizeof(float)); //A_T|(PXP)
    multiply(D_T, P, N, D, N, P, A);

    //begin QR-algorithm for A:
    float *A_current = A; //PxP; initialised to A_0
    float *E_current = initialize_identity(P); //PXP; initialised to E_0
    float *A_next = (float *) calloc(P*P, sizeof(float));
    float *E_next = (float *) calloc(P*P, sizeof(float));
    float *Q_ = (float *) calloc(P*P, sizeof(float));
    float *R_ = (float *) calloc(P*P, sizeof(float));
    float diff_norm;
    do //convergence condition for QR-algorithm
    {
        transpose(A_current, P, P, A_T);
        classicalGS(A_current, A_T, P, Q_, R_);
        multiply(R_, P, P, Q_, P, P, A_next);
        multiply(E_current, P, P, Q_, P, P, E_next);
        diff_norm = l2_norm_diagonal_diff(A_next, A_current, P);
        A_current = A_next;
        E_current = E_next;
        printf("diff_norm: %f\n", diff_norm);
    }
    while(diff_norm<TOLERANCE);
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{

    
}
