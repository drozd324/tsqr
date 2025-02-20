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
#include <mpi.h>
#include "tsqr.h"

/*
in pseudo code the serial version

alloc Q
alloc R

for 2 iter
	alloc Q_temp
	alloc R_temp

	do serial qr decomp

	mutliply Q by Q_temp and save to Q

	construct newly blocked R and save to R

alloc q_temp

do serial qr decomp -> save to final Q, r
mutliply Q by Q_temp and save to q

*/
/**
 * @brief Computes a QR factorisation of an input matrix struct a, using the TQSR 
 *        (Tall Skinny QR) method in a parallel way using 4 nodes.
 * 
 * @param a Pointer to desired matrix struct to deconstruct
 * @param q Pointer to an empty matrix struct
 * @param r Pointer to an empty matrix struct
*/
void tsqr(matrix* a, matrix* q, matrix* r){
    int proc;
    int num_procs;
    MPI_Status status;
    MPI_Init(NULL, NULL);
    MPI_Comm_proc(MPI_COMM_WORLD, &proc);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 

    int max_block_rows = 4;
    int rows_num = max_block_rows;
		
	// INIT FIRST LOCAL MATRICES FOR DECOMPOSITION
    if (proc == 0){

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

        // init Q to store blocks from local QR decomp
        block_matrix Q_temp;
        Q_temp.m = rows_num;
        Q_temp.n = rows_num;
        Q_temp.block_mat = malloc(rows_num * rows_num * sizeof(matrix)); //
        for (int i=0; i<rows_num; i++){ //row
            for (int j=0; j<rows_num; j++){ //col
                ((Q_temp.block_mat)[i*(Q_temp.n) + j]).mat = calloc((((R.block_mat)[i]).m) * ((R.block_mat)[i]).n, sizeof(double)); //
                ((Q_temp.block_mat)[i*(Q_temp.n) + j]).m   = ((R.block_mat)[i]).m;
                ((Q_temp.block_mat)[i*(Q_temp.n) + j]).n   = ((R.block_mat)[i]).n;
            }
        }
     
        // send to each processor except for proc=0
        // distribute a to local processors
        for (int i=1; i<num_procs; i++){
            MPI_Send(((R.block_mat)[i]).m, 1 , MPI_INT, 101+(i*10), MPI_COMM_WORLD);
            MPI_Send(((R.block_mat)[i]).n, 1, MPI_INT, 102+(i*10), MPI_COMM_WORLD);
            MPI_Send(((R.block_mat)[i]).mat, ((R.block_mat)[i]).m * ((R.block_mat)[i]).n, MPI_DOUBLE, 100+(i*10), MPI_COMM_WORLD);
        }
    }

    // recieve/distribute the deomposed bits of R into local memory in each processor
    matrix loc_R;
    if (proc != 0){
        MPI_Recv(loc_R.m, 1, MPI_INT, 101+(proc*10), MPI_COMM_WORLD, &status);
        MPI_Recv(loc_R.n, 1, MPI_INT, 102+(proc*10), MPI_COMM_WORLD, &status);
        loc_R.mat = malloc(loc_R.m * loc_R.n * sizeof(double));
        MPI_Recv(loc_R.mat, loc_R.m * loc_R.n, MPI_DOUBLE, 100+(proc*10), MPI_COMM_WORLD, &status);
    } else {  // on proc 1
        loc_R.m   = ((R.block_mat)[0]).m; 
        loc_R.n   = ((R.block_mat)[0]).n;
        loc_R.mat = ((R.block_mat)[0]).mat;
    }


    // main iteration loop
    for (int iter=0; iter<2; iter++){
	
		// DO LOCAL QR DECOMPOSITION
        if ((iter == 0) || (iter == 1 && (proc % 2))){
            matrix loc_Q_temp;
            loc_Q_temp.m   = loc_R.m;
            loc_Q_temp.n   = loc_R.n;
            loc_Q_temp.mat = calloc(loc_Q_temp.m * loc_Q_temp.n, sizeof(double)); //
                
            matrix loc_R_temp;
            loc_R_temp.m   = loc_R.n;
            loc_R_temp.n   = loc_R.n;
            loc_R_temp.mat = calloc(loc_R_temp.n * loc_R_temp.n, sizeof(double)); //

            qr_decomp(
                loc_R.m,
                loc_R.n,
                loc_Q_temp.mat,
                loc_R_temp.mat
            ); 
        }	

	    // SEND loc_r_temp TO "NEIGHBOUR" PROCESSOR 
	    if ( ((iter == 0) && (1 == (proc%2))) || ((iter == 1) && (proc == 2)) ){
	        MPI_Send(loc_R_temp.m  , 1                          , MPI_INT   , 301+((proc-1)*10), MPI_COMM_WORLD);
	        MPI_Send(loc_R_temp.n  , 1                          , MPI_INT   , 302+((proc-1)*10), MPI_COMM_WORLD);
	        MPI_Send(loc_R_temp.mat, loc_R_temp.m * loc_R_temp.n, MPI_DOUBLE, 300+((proc-1)*10), MPI_COMM_WORLD);
             
	    } 
		
		// RECIEVE FROM NEIGHBOUR PROCESSOR 
		if ( ((iter == 0) && (0 == (proc%2))) || ((iter == 1) && (proc == 0 )) ){
            // i can probably rewrite this to not have to bother with loc_temp_block    
			// init temporary block matrix
            block_matrix loc_temp_block;
            loc_temp_block.m = 2;
            loc_temp_block.n = 1;
            loc_temp_block.block_mat = malloc(loc_temp_block.m * loc_temp_block.n * sizeof(matrix));
			
		    // write local computed loc_R from proc
            ((temp_block.block_mat)[0]).mat = loc_R_temp.mat;
            ((temp_block.block_mat)[0]).m   = loc_R_temp.m;
            ((temp_block.block_mat)[0]).n   = loc_R_temp.n;
               
	        MPI_Recv(loc_R_temp.m  , 1                          , MPI_INT   , 301+((proc)*10), MPI_COMM_WORLD, %status);
	        MPI_Recv(loc_R_temp.n  , 1                          , MPI_INT   , 302+((proc)*10), MPI_COMM_WORLD, %status);
	        MPI_Recv(loc_R_temp.mat, loc_R_temp.m * loc_R_temp.n, MPI_DOUBLE, 300+((proc)*10), MPI_COMM_WORLD, %status);

		    // rewrite the local R matrix
            free(loc_R.mat);
		    matrix loc_R;
		    loc_R.m = ((temp_block.block_mat)[0]).m + ((temp_block.block_mat)[1]).m;
		    loc_R.n = ((temp_block.block_mat)[0]).n;
		    loc_R.mat = malloc(loc_R.m * loc_R.n * sizeof(double));
		
            comp_matrix(&loc_temp_block, &loc_R);
		    free(loc_temp_block.block_mat[1].mat);
		    free(loc_temp_block.block_mat);
		}

		if (proc != 0){
            // send loc_Q_temp
            MPI_Send(loc_Q_temp.m  , 1                          , MPI_INT   , 201+(proc*10), MPI_COMM_WORLD);
            MPI_Send(loc_Q_temp.n  , 1                          , MPI_INT   , 202+(proc*10), MPI_COMM_WORLD);                    
            MPI_Send(loc_Q_temp.mat, loc_Q_temp.m * loc_Q_temp.n, MPI_DOUBLE, 200+(proc*10), MPI_COMM_WORLD);

		} else { //if i==0
			// collect Q from each processor that has one
		    for (int j=1; j<rows_num; j++){
	            MPI_Recv(Q_temp.block_mat[j*rows_num + j].m  , 1                                                                      , MPI_INT   , 201+((j)*10), MPI_COMM_WORLD, &status);
	            MPI_Recv(Q_temp.block_mat[j*rows_num + j].n  , 1                                                                      , MPI_INT   , 202+((j)*10), MPI_COMM_WORLD, &status);
	            MPI_Recv(Q_temp.block_mat[j*rows_num + j].mat, Q_temp.block_mat[j*rows_num + j].m * Q_temp.block_mat[j*rows_num + j].n, MPI_DOUBLE, 200+((j)*10), MPI_COMM_WORLD, &status);

		    	Q_temp.block_mat[j*rows_num + j].m   = loc_Q_temp.m;
            	Q_temp.block_mat[j*rows_num + j].n   = loc_Q_temp.n;
            	Q_temp.block_mat[j*rows_num + j].mat = loc_Q_temp.mat;
		    }
		
		    // to make matrix out of block matrix
            matrix Q_temp_matrix;
            Q_temp_matrix.m = ((rows_num * 2) == max_block_rows) ? (a->m) : (rows_num * 2) * 2 * (a->n);
            Q_temp_matrix.n = rows_num * 2 * (a->n);
            Q_temp_matrix.mat = malloc(Q_temp_matrix.m * Q_temp_matrix.n * sizeof(double)); //*
            comp_matrix(&Q_temp, &Q_temp_matrix); // makes the block_matrix struct Q_temp into a matrix struct Q_temp_matrix
            
            // computing partial product of Q
            matrix mat_c;
            mat_c.mat = malloc(Q.m * Q_temp_matrix.n * sizeof(double)); //*
            mat_mul(&Q, &Q_temp_matrix, &mat_c);

            free(Q_temp_matrix.mat); //*
            
            Q.n = Q_temp_matrix.n;
            free(Q.mat); //*
    	    Q.mat = mat_c.mat;
		
	        rows_num = rows_num / 2;
		}
	    
     
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

    MPI_Finalize();
}
