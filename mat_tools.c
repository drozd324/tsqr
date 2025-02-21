/**
 * @file mat_tools.c 
 * @brief Contains functions for dealing with matices in array form and 
 *        the matrix and block_matrix structs.
 *
 * @author Patryk Drozd
 */

#include <stdio.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>
#include <string.h>
#include "mat_tools.h"

/**
 * @brief Prints a matrix from a single contigious array
 *
 * @param m Number of rows of matrix
 * @param n Number of columns of matrix
 * @param a Pointer to array of doubles representing matrix
 * @param lda Leading dimension of the matrix 
*/
void print_matrix(int m, int n, double* a, int lda) {
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++){
            printf("%2.2f ", a[i*lda + j]);
        }
        printf("\n");
    }
}

/**
 * @brief Prints each matrix within the block_matrix 
 *        struct
 *
 * @param block_a Poiter to a block_matrix struct with every 
 *        parameter filled out
*/
void print_block_matrix(block_matrix* block_a){
    for (int i=0; i<block_a->m; i++){
        for (int j=0; j<block_a->n; j++){
            print_matrix(block_a->block_mat[i*block_a->n + j].m,
                         block_a->block_mat[i*block_a->n + j].n,
                         block_a->block_mat[i*block_a->n + j].mat,
                         block_a->block_mat[i*block_a->n + j].n
            );
            printf("\n");
        }
    }
}


/**
 * @brief A function which combines two of LAPACKE's qr factorisation
 *        function to obain 
 *
 * @param m Number of rows of matrix to decompose
 * @param n Number of columns of matrix to decompose
 * @param a Pointer to array of doubles representing matrix to decompose
 * @param q Pointer to array of doubles representing matrix to put Q
 *          from QR = A into.
 * @param r Pointer to array of doubles representing matrix to put R
 *          from QR = A into.
*/
void qr_decomp(int m, int n, double* a, double* q, double* r){
    int lda = n;
    double* tau = malloc(fmin(m, n) * sizeof(double));
	
    memcpy(q, a, n*m * sizeof(double));
	
    // Double precision GEneral QR Factorisation
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, m, n, q, lda, tau); 
	
    // write to R to r
    for (int i=0; i<n; i++) {
        for (int j=i; j<n; j++){
            r[i*n + j] = q[i*n + j];
        }
    }

    int k = fmin(m, n);
    // Double precision ORthogonal Generate? QR Factorisation
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, m, n, k, q, lda, tau); 
    free(tau);
}


/**
 * @brief Condensed version of cblas dgemm so that it is easier to read 
 *        when called in code. Performs a matrix multiplication on the 
 *        matrices a and b. 
 * 
 * @param a Pointer to matrix struct a with m,n
 * @param b Pointer to matrix struct b with m,n
 * @param c Pointer to result matrix
*/
void mat_mul(matrix* a, matrix* b, matrix* c){
    int alpha = 1;
    int beta = 0;
 
    c->m = a->m;
    c->n = b->n;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        a->m, b->n, a->n, 
	alpha, 
        a->mat, a->n,
        b->mat, b->n, 
        beta,
	c->mat, c->n);
}


/**
 * @brief Takes in the lenght of an array n, the desired amount of blocks 
 *        to split it into N, the index of the desired block block_index,
 *        and gives the user start and endpoints for this block int the 
 *        original array s, e.
 * 
 * @param n Number of entries in array 
 * @param N Number of desired blocks 
 * @param block_index What would be the index if the block 
 * @param s Start of block indexed by block_index in original array
 * @param e End of block indexed by block_index in original array
*/
void decomp1d(int n, int N, int block_index, int* s, int* e){
    int remainder = n % N;
    int base = n / N;

    if (block_index < remainder){
        *s =  block_index * (base+1);
        *e = (*s) + base+1;
    } else {
        *s = (remainder * (base+1)) + ((fmax(block_index-remainder, 0)) * base);
        *e = (*s) + base;
    }
}

//void decomp_matrix(matrix* a, int M, int N, block_matrix* block_a){
//block_a->m = M;
//block_a->n = N;	
//block_a->block_matrix = malloc(block->m * block->n * sizeof(matrix*));

/**
 * @brief Decomposes an input matrix struct into a block_matrix struct
 *        with dimensions of block matrix given in the block_matrix 
 *        struct.
 * 
 * @param a Pointer to an empty matrix struct
 * @param block_a Pointer to a block_matrix struct with assigned shape, (n and m).
*/
void decomp_matrix(matrix* a, block_matrix* block_a){
    int s1;
    int e1;
    int s2;
    int e2;
    
    // for each desired block
    for (int i=0; i < (block_a->m); i++){
        decomp1d(a->m, block_a->m, i, &s1, &e1);
        for (int j=0; j < (block_a->n); j++){
            decomp1d(a->n, block_a->n , j, &s2, &e2);

            ((block_a->block_mat)[i*(block_a->n) + j]).m   = e1 - s1; 
            ((block_a->block_mat)[i*(block_a->n) + j]).n   = e2 - s2; 
            ((block_a->block_mat)[i*(block_a->n) + j]).mat = malloc((e1 - s1) * (e2 - s2) * sizeof(double));
            // gotta free this memory yourself after 
            // you're finished with the block matrix
            
            // iterate over each entry of this block and plop in corresponding 
            // values from the matrix struct a
            for (int k=0; k < (e1-s1); k++){
                for (int l=0; l < (e2-s2); l++){
                    (((block_a->block_mat)[i*(block_a->n) + j]).mat)[k*(e2-s2) + l] = (a->mat)[(s1+k)*(a->n) + s2+l]; 
                }
            }
        }
    }
}


/**
 * @brief frees allocated memory by the decomp_matrix function. 
 * 
 * @param block_a pointer to a block_matrix struct with assigned shape, (n and m).
*/
void free_block_matrix(block_matrix* block_a){
    // iterate over each entry of block matrix
    for (int i=0; i<block_a->m; i++){
        for (int j=0; j<block_a->n; j++){
            // picks out address of array within matrix struct within block matrix to free
            free(((block_a->block_mat)[i*(block_a->n) + j] ).mat);
        }
    }
    
    free(block_a->block_mat);
}


//void comp_matrix(int m, int n, double* a, int M, int N, double* a){
/**
 * @brief Reconstructs (composes) a matrix struct back from a block_block matrix struct
 *        created by decomp_matrix.
 * 
 * @param block_a Pointer to a block_matrix struct with assigned shape, (n and m).
 * @param a Pointer to an empty matrix struct. This needs to have allocated entris for the matrix
*/
void comp_matrix(block_matrix* block_a, matrix* a){
    int row_step = 0;
    int col_step = 0;
    int sub_rows = 0;
    int sub_cols = 0;
	
    // iterating over blocks
    for (int i=0; i<(block_a->m); i++){
        for (int j=0; j<(block_a->n); j++){
            sub_rows = ((block_a->block_mat)[i*(block_a->n) + j]).m;
            sub_cols = ((block_a->block_mat)[i*(block_a->n) + j]).n;
	    
            // iterating over entries in block
            for (int k=0; k<sub_rows; k++){
                for (int l=0; l<sub_cols; l++){
                    (a->mat)[(row_step + k)*(a->n) + (col_step + l)] = (((block_a->block_mat)[i*(block_a->n) + j]).mat)[k*sub_cols + l];
                }
            }

            // picks out address of array within matrix struct, within block matrix to free
            col_step += sub_cols;
        }

        row_step += sub_rows;
        col_step = 0;
    }
}
