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

    srand(123);
    for (int i=0; i<(m*n); i++){
        (a.mat)[i] = (double)rand() / (double)RAND_MAX;
    }
    
    //print_matrix(int m, int n, double* a, int lda);
    printf("printing input matrix a\n");
    print_matrix(a.m, a.n, a.mat, a.n);

    //decomp1d(int n, int N, int block_index, int s, int e);
 
    block_matrix block_a;
    block_a.m = 3;
    block_a.n = 1;
    block_a.block_mat = malloc((block_a.m) * (block_a.n) * sizeof(matrix));
 
    //decomp_matrix(matrix* a, block_matrix* block_a);
    printf("decomp_mat start\n");
    decomp_matrix(&a, &block_a);
    printf("decomp_mat end\n");

    printf("printing block_a\n");
    for (int i=0; i<block_a.m; i++){
        for (int j=0; j<block_a.n; j++){
            print_matrix(block_a.block_mat[i*block_a.n + j].m,
                         block_a.block_mat[i*block_a.n + j].n,
                         block_a.block_mat[i*block_a.n + j].mat,
                         block_a.block_mat[i*block_a.n + j].n
            );
            printf("\n");
        }
    }

    //comp_matrix(block_matrix* block_a, matrix* a);

    matrix comp_a;
    comp_a.mat = malloc(m * n * sizeof(double));
    comp_a.m = m;
    comp_a.n = n;
    printf("comp_mat start\n");
    comp_matrix(&block_a, &comp_a);
    printf("comp_mat end\n");

    printf("printing input matrix comp_a\n");
    print_matrix(comp_a.m, comp_a.n, comp_a.mat, comp_a.n);
    
    //free_block_matrix(block_matrix* block_a);


    return 0; 
}

