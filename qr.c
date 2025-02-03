#include <stdio.h>
#include <lapacke.h>

void print_matrix(const char* desc, int m, int n, double* a, int lda) {
    printf("\n %s\n", desc);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++){
            printf("%8.4f ", a[i * lda + j]);
		}
        printf("\n");
    }
}

int main(){
	int m = 5, n = 3;
	double* tau = malloc(min(m,n) * sizeof(double));
	double* a = malloc(m * n * sizeof(double));
    
    int lda = m;
 
	LAPACK_dgemm(); // Double precision GEneral Matrix Matrix
    LAPACKE_dgeqrf(m, n, *a, lda, tau, work, lwork, info); // Double precision GEneral QR Factorisation
    LAPACKE_dorgqr(m, n, *a, lda, tau, work, lwork, info); // Double precision ORthogonal Generate? QR Factorisation

    return info; 
}
