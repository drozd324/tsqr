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

    // initialises q matrix to be an identity matrix
    for (int i=0; i<(a->m); i++){
        (q->mat)[i*(a->n) + i] = 1;
    }
    
    // initialise block matrix R for a to be decomposed into 
    block_matrix* R;
    R->m = rows_num;
    R->n = 1;
    decomp_matrix(a, R);
	
    // main iteration loop
    for (int iter=0; iter<3; iter++){
        
        // init temp matrices for computations in current iteration		
        // init Q		
        block_matrix* Q_temp;
        Q_temp->m = rows_num;
        Q_temp->n = rows_num;
        Q_temp->block_mat = malloc(rows_num * rows_num * sizeof(matrix));
        
        // init R		
        block_matrix* R_temp;
        R_temp->m = rows_num;
        R_temp->n = 1;
        R_temp->block_mat = malloc(rows_num * sizeof(matrix));
        
        // allocate memory for Q_temp
        for (int i=0; i<rows_num; i++){
            for (int j=0; j<rows_num; j++){
                // (array of matrix structs)[of matrix struct]
                // (block matrix)[index of matrix]->(array of double representing matrix)  =  calloc size m_i*n array of doubles representing matrix
                ((Q_temp->block_mat)[i*(Q_temp->n) + j]).mat = calloc(((R->block_mat)->m) * (a->n), sizeof(double*));
                ((Q_temp->block_mat)[i*(Q_temp->n) + j]).m   = (R->block_mat)->m;
                ((Q_temp->block_mat)[i*(Q_temp->n) + j]).n   = a->n;
            }
        } 
        
        // allocate memory for R_temp
        for (int i=0; i<rows_num; i++){
            ((R_temp->block_mat)[i]).mat = calloc((a->n) * (a->n), sizeof(double));
            ((R_temp->block_mat)[i]).m   = (R->block_mat)->m;
            ((R_temp->block_mat)[i]).n   = a->n;
        }
        
        // perform qr decomp for each block		
        for (int i=0; i<rows_num; i++){
            // index the block matrix and pick out the array of doubles representing matrix 
            //        m  n  (mat to decomp)
            qr_decomp(
                ((R->block_mat)[i]).m,                         // m
                ((R->block_mat)[i]).n,                         // n 
                ((R->block_mat)[i]).mat,                       // matrix to decompose
                ((Q_temp->block_mat)[i*(Q_temp->n) + i]).mat,  // pick out  (i,i) entry in Q_temp which would be a matrix
                ((R_temp->block_mat)[i]).mat                   // pick out (i,0) entry in block mat which is again a matrix
            );                                                                    
        } 
        
        // matrix needed for computation
        matrix* Q_temp_matrix;
        Q_temp_matrix->m = (rows_num == max_block_rows) ? (a->m) : rows_num * 2 * (a->n);
        Q_temp_matrix->n = rows_num * (a->n);
        comp_matrix(Q_temp, Q_temp_matrix); // makes the block_matrix struct Q_temp into a matrix struct Q_temp_matrix
	
        matrix* mat_c;
        mat_c->mat = calloc(q->m * Q_temp->n, sizeof(double));
        mat_mul(q, Q_temp_matrix, mat_c);
        memcpy(q->mat, mat_c->mat, (mat_c->m) * (mat_c->n) * sizeof(double));
    
        rows_num = rows_num / 2;
        free_block_matrix(R);
        R->m = rows_num;
         
        for (int i=0; i<rows_num; i++){
            block_matrix* temp_block;
            temp_block->m = 2;
            temp_block->n = 1;
            ((temp_block->block_mat)[0]).mat = malloc( (((R_temp->block_mat)[i*2    ]).m) * (R->n) * sizeof(double));
            ((temp_block->block_mat)[1]).mat = malloc( (((R_temp->block_mat)[i*2 + 1]).m) * (R->n) * sizeof(double));
           
            //comp_matrix((R->block_mat)[i].mat, temp_mat);
            comp_matrix(temp_block, &((R->block_mat)[i]));
        }

    }
    
    matrix* q_temp;
    q_temp->m = ((R->block_mat)[0]).m;
    q_temp->n = ((R->block_mat)[0]).n;
    q_temp->mat = malloc( (q_temp->m) * (q_temp->n) * sizeof(double));
    
    qr_decomp(
        ((R->block_mat)[0]).m, 
        ((R->block_mat)[0]).n,
        ((R->block_mat)[0]).mat,
        q_temp->mat,
        r->mat
    );

    matrix* mat_c;
    mat_c->mat = malloc(q->m * q->n * sizeof(double));
    
    mat_mul(q, q_temp, mat_c);
    memcpy(q->mat, mat_c->mat, (mat_c->m) * (mat_c->n) * sizeof(double)); 
}
