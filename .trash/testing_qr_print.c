//row major ideology
// electric fence for debugging!

#include "mat_tools.h"
#include <cblas.h>
#include <stdlib.h>

int main(){
    int m = 9;
    int n = 5;

    matrix a;
    a.mat = malloc(m * n * sizeof(double));
    a.m = m;
    a.n = n;

    matrix q;
    q.mat = malloc(m * n * sizeof(double));
    q.m = m;
    q.n = n;

    matrix r;
    r.mat = malloc(n * n * sizeof(double));
    r.m = n;
    r.n = n;

    matrix qr;
    qr.mat = malloc(m * n * sizeof(double));
    qr.m = m;
    qr.n = n;

    srand(1234);
    for (int i=0; i<(m*n); i++){
        (a.mat)[i] = (double)rand() / (double)RAND_MAX;
    }
   
    
    //print_matrix(int m, int n, double* a, int lda);
    printf("printing input matrix a\n");
    print_matrix(a.m, a.n, a.mat, a.n);

    qr_decomp(m, n, a.mat, q.mat, r.mat);
    
    mat_mul(&q, &r, &qr);
		
    printf("printing input matrix qr\n");
    print_matrix(qr.m, qr.n, qr.mat, qr.n);

    //decomp1d(int n, int N, int block_index, int s, int e);
 
    block_matrix block_a;
    block_a.m = 3;
    block_a.n = 2;
    block_a.block_mat = malloc((block_a.m) * (block_a.n) * sizeof(matrix*));
 
    //decomp_matrix(matrix* a, block_matrix* block_a);
    printf("decomp_mat start");
    decomp_matrix(&a, &block_a);
    printf("decomp_mat end");

    //comp_matrix(block_matrix* block_a, matrix* a);

    matrix comp_a;
    comp_a.mat = malloc(m * n * sizeof(double));
    comp_a.m = m;
    comp_a.n = n;
    printf("comp_mat start");
    comp_matrix(&block_a, &comp_a);
    printf("comp_mat end");

    printf("printing input matrix comp_a\n");
    print_matrix(qr.m, comp_a.n, comp_a.mat, comp_a.n);
    
    //free_block_matrix(block_matrix* block_a);


    return 0; 
}

