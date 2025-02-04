#include <stdio.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>

void print_matrix(int m, int n, double* a, int lda) {
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++){
            printf("%2.2f ", a[i*lda + j]);
		}
        printf("\n");
    }
}

int main(){
	srand(1234);
	
	int m = 4;
	int n = 3;
	double* a = malloc(m * n * sizeof(double));
	for (int i=0; i<(m*n); i++){
		//printf("%d\n", i);
		a[i] = (double)rand() / (double)RAND_MAX;
	}
	
	double* r = calloc(n * n, sizeof(double));

    int     lda   = n;
	double* tau   = malloc(fmin(m, n) * sizeof(double));
	int 	info;

	printf("printing a before \n"); 
	print_matrix(m, n, a, lda);
	printf("printing r before \n"); 
	print_matrix(n, n, r, lda);

	// Double precision GEneral QR Factorisation
	printf("LAPACK_ROW_MAJOR = %d", LAPACK_ROW_MAJOR);
	printf("LAPACK_COL_MAJOR = %d", LAPACK_COL_MAJOR);
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, m, n, a, lda, tau); 
	
	// write to R to r
    for (int i=0; i<n; i++) {
        for (int j=i; j<n; j++){
			//printf("%d\n", i*n + j);
            r[i*n + j] = a[i*n + j];
		}
	}

	printf("\n");
	printf("printing a after r \n"); 
	print_matrix(m, n, a, lda);
	printf("printing r after r \n"); 
	print_matrix(n, n, r, lda);

	int k = fmin(m, n);
	// Double precision ORthogonal Generate? QR Factorisation
   	LAPACKE_dorgqr(LAPACK_ROW_MAJOR, m, n, k, a, lda, tau); 

	printf("\n");
	printf("printing a after q \n"); 
	print_matrix(m, n, a, lda);
	printf("printing r after q \n"); 
	print_matrix(n, n, r, lda);
	
	// Perform matrix multiplication C = alpha * A * B + beta * C
	double alpha = 1;
	double beta = 0;
	
	double* c = calloc(m * n, sizeof(double));
	// THIS NO WORKIE
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, a[0][0], lda, r[0][0], lda, beta, &c[0][0], lda);

	printf("\n");
	printf("printing c, oringinal a\n"); 
	print_matrix(m, n, c, lda);

    return 0; 
}
