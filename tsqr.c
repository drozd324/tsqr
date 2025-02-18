/**
 * @file tsqr.c
 * @brief TSQR function for qr factorisation
 *
 * @author Patryk Drozd
 */

#include <stdio.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>
#include <string.h>
#include "tsqr.h"

////////////////////////////////////////////////////////////////////////////////////////////////
////CHECK FOR MEMORY LEAKS//////////CHECK FOR MEMORY LEAKS//////////CHECK FOR MEMORY LEAKS//////
////////////////////////////////////////////////////////////////////////////////////////////////

//UPDATE SHAPES OF MATRICES AND BLOCK MATRICES IN THE CODE//

/**
 * @brief Computes a QR factorisation of an input matrix struct a, using the TQSR 
 *        (Tall Skinny QR) method in a parallel way using 4 nodes.
 * 
 * @param a Pointer to desired matrix struct to deconstruct
 * @param q Pointer to an empty matrix struct
 * @param r Pointer to an empty matrix struct
*/
void tsqr(matrix* a, matrix* q, matrix* r){
    int max_block_rows = 4;
    int rows_num = max_block_rows;
    //int m = a->m;
    //int n = a->n;

    // initialises Q matrix to be an identity matrix
    matrix Q;
    Q.m = q->m;
    Q.n = q->m;
    Q.mat = calloc(Q.m * Q.m, sizeof(double)); //
    for (int i=0; i<(Q.m); i++){
        (Q.mat)[i*(Q.m) + i] = 1;
    }
    
    // initialise block matrix R for a to be decomposed into 
    block_matrix R;
    R.m = rows_num;
    R.n = 1;
    R.block_mat = malloc(R.m * R.n * sizeof(matrix));
    decomp_matrix(a, &R); //

    // main iteration loop
    for (int iter=0; iter<2; iter++){
	printf("PASS %d\n", iter);
        // init temp matrices for computations in current iteration		
        // init Q		
        block_matrix Q_temp;
        Q_temp.m = rows_num;
        Q_temp.n = rows_num;
        Q_temp.block_mat = malloc(rows_num * rows_num * sizeof(matrix)); //
        
        // init R	
        block_matrix R_temp;
        R_temp.m = rows_num;
        R_temp.n = 1;
        R_temp.block_mat = malloc(rows_num * sizeof(matrix)); //
        
        // allocate memory for Q_temp
        for (int i=0; i<rows_num; i++){ //row
            for (int j=0; j<rows_num; j++){ //col
                // (array of matrix structs)[of matrix struct]
                // (block matrix)[index of matrix]->(array of double representing matrix)  =  calloc size m_i*n array of doubles representing matrix
                ((Q_temp.block_mat)[i*(Q_temp.n) + j]).mat = calloc((((R.block_mat)[i]).m) * ((R.block_mat)[i]).n, sizeof(double)); //
                ((Q_temp.block_mat)[i*(Q_temp.n) + j]).m   = ((R.block_mat)[i]).m;
                ((Q_temp.block_mat)[i*(Q_temp.n) + j]).n   = ((R.block_mat)[i]).n;
            }
        } 
        
        // allocate memory for R_temp
        for (int i=0; i<rows_num; i++){
            ((R_temp.block_mat)[i]).mat = calloc((R.block_mat[i]).m * (R.block_mat[i]).n, sizeof(double)); //
            ((R_temp.block_mat)[i]).m   = (R.block_mat[i]).m;
            ((R_temp.block_mat)[i]).n   = (R.block_mat[i]).n;
        }
    	
        // perform qr decomp for each block		
	// add check for m>n so that it doesnt shit itself
        for (int i=0; i<rows_num; i++){
            // index the block matrix and pick out the array of doubles representing matrix 
            //        m  n  (mat to decomp)
            //printf("m_i %d\n", ((R.block_mat)[i]).m);
            //printf("n_i %d\n", ((R.block_mat)[i]).n);
	    //print_matrix(((R.block_mat)[i]).m, ((R.block_mat)[i]).n, ((R.block_mat)[i]).mat, ((R.block_mat)[i]).n);
            qr_decomp(
                ((R.block_mat)[i]).m,                         // m
                ((R.block_mat)[i]).n,                         // n 
                ((R.block_mat)[i]).mat,                       // matrix to decompose
                ((Q_temp.block_mat)[i*(Q_temp.n) + i]).mat,   // pick out  (i,i) entry in Q_temp which would be a matrix
                ((R_temp.block_mat)[i]).mat                   // pick out (i,0) entry in block mat which is again a matrix
            );
	    //print_matrix(((R.block_mat)[i]).m, ((R.block_mat)[i]).n, ((R.block_mat)[i]).mat                    , ((R.block_mat)[i]).n);
	    //print_matrix(((R.block_mat)[i]).m, ((R.block_mat)[i]).n, ((Q_temp.block_mat)[i*(Q_temp.n) + i]).mat, ((R.block_mat)[i]).n);
	    //print_matrix(((R.block_mat)[i]).m, ((R.block_mat)[i]).n, ((R_temp.block_mat)[i]).mat               , ((R.block_mat)[i]).n);
        
        } 

        // matrix needed for computation
        matrix Q_temp_matrix;
        Q_temp_matrix.m = (rows_num == max_block_rows) ? (a->m) : rows_num * 2 * (a->n);
        Q_temp_matrix.n = rows_num * (a->n);
        Q_temp_matrix.mat = malloc(Q_temp_matrix.m * Q_temp_matrix.n * sizeof(double)); //*
	
	printf("Q.m = %d\n", Q.m);
	printf("Q.n = %d\n", Q.n);
	printf("Q_temp.m = %d\n", Q_temp.m);
	printf("Q_temp.n = %d\n", Q_temp.n);
	print_matrix(Q.m, Q.n, Q.mat, Q.m);
	
	//print_block_matrix(&Q_temp);
      
        comp_matrix(&Q_temp, &Q_temp_matrix); // makes the block_matrix struct Q_temp into a matrix struct Q_temp_matrix
	
	// computing partial product of Q
        matrix mat_c;
        mat_c.mat = malloc(Q.m * Q_temp_matrix.n * sizeof(double)); //*
	//printf("Q.m, Q.n = %d, %d\n", Q.m, Q.n);
	printf("Q_temp_matrix.m, Q_temp_matrix.n = %d, %d\n", Q_temp_matrix.m, Q_temp_matrix.n);
        mat_mul(&Q, &Q_temp_matrix, &mat_c);
        free(Q_temp_matrix.mat); //*
        
        //Q.m = q->m; stays the same
        Q.n = Q_temp_matrix.n;
	free(Q.mat); //*
        Q.mat = malloc(Q.m * Q.n * sizeof(double)); //

	memcpy(Q.mat, mat_c.mat, Q.m * Q.n * sizeof(double));
        free(mat_c.mat); //*
    
        rows_num = rows_num / 2;
	//printf("R before\n");
	//print_block_matrix(&R);
        free_block_matrix(&R);

        R.m = rows_num;
        R.n = 1;
	//printf("R.m, R.n = %d, %d\n", R.m, R.n);
    	R.block_mat = malloc(R.m * R.n * sizeof(matrix));
         
        for (int i=0; i<rows_num; i++){
            printf("i=%d\n", i);
            block_matrix temp_block;
            temp_block.m = 2;
            temp_block.n = 1;
	    temp_block.block_mat = malloc(temp_block.m * temp_block.n * sizeof(matrix));
		
	    printf("printing the two block that are to be putinto temp_block\n");
	    print_matrix(((R_temp.block_mat)[i*2    ]).m, ((R_temp.block_mat)[i*2    ]).n, ((R_temp.block_mat)[i*2    ]).mat, ((R_temp.block_mat)[i*2    ]).n);
	    printf("\n");
	    print_matrix(((R_temp.block_mat)[i*2 + 1]).m, ((R_temp.block_mat)[i*2 + 1]).n, ((R_temp.block_mat)[i*2 + 1]).mat, ((R_temp.block_mat)[i*2 + 1]).n);

            ((temp_block.block_mat)[0]).mat = ((R_temp.block_mat)[i*2    ]).mat;
            ((temp_block.block_mat)[1]).mat = ((R_temp.block_mat)[i*2 + 1]).mat;
            ((temp_block.block_mat)[0]).m = ((R_temp.block_mat)[i*2    ]).m;
            ((temp_block.block_mat)[1]).m = ((R_temp.block_mat)[i*2 + 1]).m;
            ((temp_block.block_mat)[0]).n = ((R_temp.block_mat)[i*2    ]).n;
            ((temp_block.block_mat)[1]).n = ((R_temp.block_mat)[i*2 + 1]).n;

            printf("printing temp block for i=%d\n", i);
	    print_block_matrix(&temp_block);

            ((R.block_mat)[i]).m = ((R_temp.block_mat)[i*2]).m + ((R_temp.block_mat)[i*2 + 1]).m;
            ((R.block_mat)[i]).n = a->n;
            R.block_mat[i].mat = malloc( (R.block_mat[i]).m * (R.block_mat[i]).n * sizeof(double));

            comp_matrix(&temp_block, &((R.block_mat)[i])); 
        }
	    
	printf("R after\n");
	print_block_matrix(&R);

        free_block_matrix(&Q_temp);
        free_block_matrix(&R_temp);

    }
    
    matrix q_temp;
    q_temp.m = ((R.block_mat)[0]).m;
    q_temp.n = ((R.block_mat)[0]).n;
    q_temp.mat = malloc( (q_temp.m) * (q_temp.n) * sizeof(double));
    
    qr_decomp(
        ((R.block_mat)[0]).m, 
        ((R.block_mat)[0]).n,
        ((R.block_mat)[0]).mat,
        q_temp.mat,
        r->mat
    );

    mat_mul(&Q, &q_temp, q);
    free(q_temp.mat);
}
