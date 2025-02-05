#include <stdio.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>

//row major ideology

struct _Matrix{
	double* matrix;
	int m;
	int n;
}matrix;

struct _block_Matrix{
	matrix* block_matrix;
	int m;
	int n;
}block_matrix;

void print_matrix(int m, int n, double* a, int lda); 
void qr_decomp(int m, int n, double* a, double* r);
void decomp1d(n, N, block_index);

int main(){
	int m = 4;
	int n = 3;

	double* a = malloc(m * n * sizeof(double));
	srand(1234);
	for (int i=0; i<(m*n); i++){
		a[i] = (double)rand() / (double)RAND_MAX;
	}
	
	double* r = calloc(n * n, sizeof(double));
	
	qr_decomp(m, n, a, r);

		
    return 0; 
}


void print_matrix(int m, int n, double* a, int lda) {
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++){
            printf("%2.2f ", a[i*lda + j]);
		}
        printf("\n");
    }
}

void qr_decomp(int m, int n, double* a, double* r){
    int     lda = n;
	double* tau = malloc(fmin(m, n) * sizeof(double));
	
	// Double precision GEneral QR Factorisation
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, m, n, a, lda, tau); 
	
	// write to R to r
    for (int i=0; i<n; i++) {
        for (int j=i; j<n; j++){
            r[i*n + j] = a[i*n + j];
		}
	}

	int k = fmin(m, n);
	// Double precision ORthogonal Generate? QR Factorisation
   	LAPACKE_dorgqr(LAPACK_ROW_MAJOR, m, n, k, a, lda, tau); 
}


void decomp1d(int n, int N, int block_index, int s, int e){
    int remainder = n % N;
    int base = n / N;

    if (block_index < remainder){
        s =  block_index * (base+1);
        e = s + base+1;
    } else {
        s = (remainder * (base+1)) + ((fmax(block_index-remainder, 0)) * base);
        e = s + base;
    }
}

void decomp_matrix(matrix* a, int M, int N, block_matrix* block_a){
	int m = a.m;
	int n = a.n;	
	double* a = a.matrix;

	block_a.m = M;
	block_a.n = N;

	int s1;
	int e1;
	int s2;
	int e2;

	for (int i=0; i<M; i++){
        decomp1d(a.m, M, i, s1, e1);
		for (int j=0; j<N; j++){
            decomp1d(a.n, N, j, s2, e2);
           
			matrix block;
			block.m = e1 - s1;
			block.n = e2 - s2;
			
			block.matrix = malloc(block.m * block.n * sizeof(double)); 
			// gotta free this memory yourself after 
			// you're finished with the block matrix
			for (int k=0; k < block.m; k++){
				for (int l=0; l < block.n; l++){
                    block.matrix[k*block.n + l] = a.matrix[(s1+k)*a.n + s2+l];
				}
			}	

			block_a.block_matrix[i*n + j] = &block

		}
	}
}



//void comp_matrix(int m, int n, double* a, int M, int N, double* a){
void comp_matrix(block_matrix block_a, matrix a){
    int row_step = 0;
    int col_step = 0;
    int sub_rows = 0;
    int sub_cols = 0;
    
	for (int i=0; i<M; i++){
		for (int j=0; j<N; j++){
            
			sub_rows = block_a.block_matrix[i*N + j].m;
			sub_cols = block_a.block_matrix[i*N + j].n;

        	for (int k=0; k<sub_rows; k++){
        		for (int l=0; l<sub_rows; l++){
                    a.matrix[(row_step + k)*n + (col_step + l)] = block[k*sub_cols + l];
                    
            col_step += sub_cols;
        row_step += sub_rows;
        col_step = 0;
	}
}

