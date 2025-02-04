#include <stdio.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>

//row major ideology

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


void decomp1d(int n, int N, int block_index, int* s, int* e){
    int remainder = n % N;
    int base = n / N;

    if (block_index < remainder){
        *s =  block_index * (base+1);
        *e = *s + base+1;
    } else {
        *s = (remainder * (base+1)) + ((fmax(block_index-remainder, 0)) * base);
        *e = *s + base;
    }
}

//double* block_a = malloc(M * N * sizeof(double*));
void decomp_matrix(int m, int n, double* a, int M, int N, double* block_a){
    
	for (int i=0; i<M; i++){
		int* s1;
		int* e1;
        decomp1d(m, M, i, s1, e1);

		for (int j=0; j<N; j++){
			int* s2;
			int* e2;
            decomp1d(n, N, j, s2, e2);
            
			&block_A[i*n + j] = malloc((e1 - s1)*(e2 - s2) * sizof(double*); 
			for (int a=0; a<(*e1 - *s1); a++){
				for (int b=0; b<(*e2 - *s2); b++){
                    block_A[i][j][a][b] = A[s1 + a][s2 + b];
				}		
			}
		}
	}
}

