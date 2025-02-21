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

/**
 * @brief Computes a QR factorisation of an input matrix struct a, using the TQSR 
 *        (Tall Skinny QR) method in parallel using 4 nodes.
 * 
 * @param a Pointer to desired matrix struct to deconstruct
 * @param q Pointer to an empty matrix struct
 * @param r Pointer to an empty matrix struct
*/
void tsqr(matrix* a, matrix* q, matrix* r){
    int proc;
    int num_procs;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 

    int max_block_rows = 4;
    int rows_num = max_block_rows;
	
    matrix Q;
    block_matrix R;
    block_matrix Q_temp;

	// INIT FIRST LOCAL MATRICES FOR DECOMPOSITION
    if (proc == 0){

        // initialises Q matrix to be an identity matrix
        //matrix Q;
        Q.m = q->m;
        Q.n = q->m;
        Q.mat = calloc(Q.m * Q.m, sizeof(double)); //
        for (int i=0; i<(Q.m); i++){
            (Q.mat)[i*(Q.m) + i] = 1;
        }
        
        // initialise block matrix R for a to be decomposed into 
        //block_matrix R;
        R.m = rows_num;
        R.n = 1;
        R.block_mat = malloc(R.m * R.n * sizeof(matrix));
        decomp_matrix(a, &R); //

        // SEND TO EACH PROCESSOR EXCEPT FOR PROC=0
        for (int i=1; i<num_procs; i++){
			//fprintf(stderr, "Send to i*10=%d from 0\n", i*10);
            MPI_Send(&(((R.block_mat)[i]).m) , 1                                          , MPI_INT   , i, 101+(i*10), MPI_COMM_WORLD);
            MPI_Send(&(((R.block_mat)[i]).n) , 1                                          , MPI_INT   , i, 102+(i*10), MPI_COMM_WORLD);
            MPI_Send(  ((R.block_mat)[i]).mat, ((R.block_mat)[i]).m * ((R.block_mat)[i]).n, MPI_DOUBLE, i, 100+(i*10), MPI_COMM_WORLD);
        }
    }

    // RECIEVE/WIRTE THE DEOMPOSED BITS OF R INTO LOCAL MEMORY IN EACH PROCESSOR
    matrix loc_R;
    if (proc != 0){
		//fprintf(stderr, "Recv to proc*10=%d from 0\n", proc*10);
        MPI_Recv(&(loc_R.m) , 1                , MPI_INT   , 0, 101+(proc*10), MPI_COMM_WORLD, &status);
        MPI_Recv(&(loc_R.n) , 1                , MPI_INT   , 0, 102+(proc*10), MPI_COMM_WORLD, &status);
        loc_R.mat = malloc(loc_R.m * loc_R.n * sizeof(double));
        MPI_Recv(  loc_R.mat, loc_R.m * loc_R.n, MPI_DOUBLE, 0, 100+(proc*10), MPI_COMM_WORLD, &status);	

    } else {  // on proc 0
        loc_R.m   = ((R.block_mat)[0]).m; 
        loc_R.n   = ((R.block_mat)[0]).n;
        loc_R.mat = ((R.block_mat)[0]).mat;
    }

    // main iteration loop
    matrix loc_Q_temp;
    matrix loc_R_temp;
    for (int iter=0; iter<2; iter++){

		// DO LOCAL QR DECOMPOSITION
        if ((iter == 0) || ((iter == 1) && (0 == (proc % 2))) ){
			//fprintf(stderr, "ITER=%d on proc=%d\n", iter, proc);

            //matrix loc_Q_temp;
            loc_Q_temp.m   = loc_R.m;
            loc_Q_temp.n   = loc_R.n;
            loc_Q_temp.mat = calloc(loc_Q_temp.m * loc_Q_temp.n, sizeof(double)); //
                
            //matrix loc_R_temp;
            loc_R_temp.m   = loc_R.n;
            loc_R_temp.n   = loc_R.n;
            loc_R_temp.mat = calloc(loc_R_temp.n * loc_R_temp.n, sizeof(double)); //

			//fprintf(stderr, "QR DECOMP ITER=%d, proc = %d \n", iter, proc);
            qr_decomp(
                loc_R.m,
                loc_R.n,
                loc_R.mat,
                loc_Q_temp.mat,
                loc_R_temp.mat
            ); 
        }	

	    // SEND loc_r_temp TO "NEIGHBOUR" PROCESSOR
	    if ( ((iter == 0) && (1 == (proc%2))) || ((iter == 1) && (proc == 2)) ){
			//fprintf(stderr, "R ITER=%d, Send to neighbour=%d from %d,   tag =  %d \n", iter, proc-(1+iter), proc, 301+((proc-(1+iter))*10));
			//fprintf(stderr, "AAAAAAA (loc_R_temp.m) = %d\n", (loc_R_temp.m));
	        MPI_Send(&(loc_R_temp.m)  , 1                          , MPI_INT   , proc-(1+iter), 301+((proc-(1+iter))*10), MPI_COMM_WORLD);
	        MPI_Send(&(loc_R_temp.n)  , 1                          , MPI_INT   , proc-(1+iter), 302+((proc-(1+iter))*10), MPI_COMM_WORLD);
	        MPI_Send(  loc_R_temp.mat , loc_R_temp.m * loc_R_temp.n, MPI_DOUBLE, proc-(1+iter), 300+((proc-(1+iter))*10), MPI_COMM_WORLD);
            																     // to proc 0 ,      300
	    } 
		
		// RECIEVE FROM "NEIGHBOUR" PROCESSOR 
		if ( ((iter == 0) && (0 == (proc%2))) || ((iter == 1) && (proc == 0)) ){
            // i can probably rewrite this to not have to bother with loc_temp_block    
			// init temporary block matrix
            block_matrix loc_temp_block;
            loc_temp_block.m = 2;
            loc_temp_block.n = 1;
            loc_temp_block.block_mat = malloc(loc_temp_block.m * loc_temp_block.n * sizeof(matrix));
			
		    // write local computed loc_R from proc
            ((loc_temp_block.block_mat)[0]).mat = loc_R_temp.mat;
            ((loc_temp_block.block_mat)[0]).m   = loc_R_temp.m;
            ((loc_temp_block.block_mat)[0]).n   = loc_R_temp.n;
               
			//
			//fprintf(stderr, "BEFORE ((loc_temp_block.block_mat)[1]).m = %d\n", ((loc_temp_block.block_mat)[1]).m);
			// this shit above baaaaaaaaaaaaad
			//fprintf(stderr, "R ITER=%d Recv on neighbour=%d from %d, tag = %d \n", iter, proc, proc + 1+iter, 301+((proc + 2*iter)*10));
	        MPI_Recv(&(((loc_temp_block.block_mat)[1]).m) , 1                                                                    , MPI_INT   , proc + 1+iter, 301+((proc*10)*(1-iter)), MPI_COMM_WORLD, &status);
	        MPI_Recv(&(((loc_temp_block.block_mat)[1]).n) , 1                                                                    , MPI_INT   , proc + 1+iter, 302+((proc*10)*(1-iter)), MPI_COMM_WORLD, &status);
			((loc_temp_block.block_mat)[1]).mat =    malloc(((loc_temp_block.block_mat)[1]).m * ((loc_temp_block.block_mat)[1]).n * sizeof(double));
	        MPI_Recv(  ((loc_temp_block.block_mat)[1]).mat, ((loc_temp_block.block_mat)[1]).m * ((loc_temp_block.block_mat)[1]).n, MPI_DOUBLE, proc + 1+iter, 300+((proc*10)*(1-iter)), MPI_COMM_WORLD, &status);

		    // rewrite the local R matrix
            free(loc_R.mat);
		    loc_R.m = ((loc_temp_block.block_mat)[0]).m + ((loc_temp_block.block_mat)[1]).m;
		    //loc_R.n = ((loc_temp_block.block_mat)[0]).n;
		    loc_R.mat = malloc(loc_R.m * loc_R.n * sizeof(double));
		
            comp_matrix(&loc_temp_block, &loc_R);
		    free(loc_temp_block.block_mat[1].mat);
		    free(loc_temp_block.block_mat);
		}
		
		// MAKING PART OF Q MATRIX BY MULTIPLYING IN WHAT WE HAVE SO FAR
		if ( ((iter == 0) && (proc != 0)) || ((iter == 1) && (proc == 2)) ){
            // send loc_Q_temp
			//fprintf(stderr, "Q ITER=%d, Send to 0 from %d,   tag =  %d \n", iter, proc, 201+(proc*10));
            MPI_Send(&(loc_Q_temp.m)  , 1                          , MPI_INT   , 0, 201+(proc*10), MPI_COMM_WORLD);
            MPI_Send(&(loc_Q_temp.n)  , 1                          , MPI_INT   , 0, 202+(proc*10), MPI_COMM_WORLD);
            MPI_Send(  loc_Q_temp.mat , loc_Q_temp.m * loc_Q_temp.n, MPI_DOUBLE, 0, 200+(proc*10), MPI_COMM_WORLD);
		} 

		if (proc == 0){
			
            // init Q to store blocks from local QR decomp
            Q_temp.m = rows_num;
            Q_temp.n = rows_num;
            Q_temp.block_mat = malloc(rows_num * rows_num * sizeof(matrix)); //

			// collect Q  from each processor that has one which becomes the diagonal part of Q_temp
			if (iter == 0){
			    for (int j=1; j<rows_num; j++){
		            MPI_Recv(&(Q_temp.block_mat[j*rows_num + j].m) , 1                                                                      , MPI_INT   , j, 201+(j*10), MPI_COMM_WORLD, &status);
		            MPI_Recv(&(Q_temp.block_mat[j*rows_num + j].n) , 1                                                                      , MPI_INT   , j, 202+(j*10), MPI_COMM_WORLD, &status);
			        ((Q_temp.block_mat)[j*(Q_temp.n) + j]).mat = malloc((Q_temp.block_mat[j*rows_num + j].m) * (Q_temp.block_mat[j*rows_num + j].n) *  sizeof(double)); 
		            MPI_Recv(  Q_temp.block_mat[j*rows_num + j].mat, Q_temp.block_mat[j*rows_num + j].m * Q_temp.block_mat[j*rows_num + j].n, MPI_DOUBLE, j, 200+(j*10), MPI_COMM_WORLD, &status);
					
		    	}
			} else { //if iter == 1
				int j=2;
		        MPI_Recv(&(Q_temp.block_mat[1*rows_num + 1].m) , 1                                                                      , MPI_INT   , j, 201+(j*10), MPI_COMM_WORLD, &status);
		        MPI_Recv(&(Q_temp.block_mat[1*rows_num + 1].n) , 1                                                                      , MPI_INT   , j, 202+(j*10), MPI_COMM_WORLD, &status);
			    ((Q_temp.block_mat)[1*(Q_temp.n) + 1]).mat = malloc((Q_temp.block_mat[1*rows_num + 1].m) * (Q_temp.block_mat[1*rows_num + 1].n) *  sizeof(double)); 
		        MPI_Recv(  Q_temp.block_mat[1*rows_num + 1].mat, Q_temp.block_mat[1*rows_num + 1].m * Q_temp.block_mat[1*rows_num + 1].n, MPI_DOUBLE, j, 200+(j*10), MPI_COMM_WORLD, &status);
			}
		
			Q_temp.block_mat[0].m   = loc_Q_temp.m;
            Q_temp.block_mat[0].n   = loc_Q_temp.n;
            Q_temp.block_mat[0].mat = loc_Q_temp.mat;

			// put zeros into the rest of Q_temp	
            for (int i=0; i<rows_num; i++){ //row
                for (int j=0; j<rows_num; j++){ //col	
					if (i != j){
        	           	((Q_temp.block_mat)[i*(Q_temp.n) + j]).mat = calloc((((Q_temp.block_mat)[i*rows_num + i]).m) * (a->n) , sizeof(double)); //
            	       	((Q_temp.block_mat)[i*(Q_temp.n) + j]).m   = (((Q_temp.block_mat)[i*rows_num + i]).m);
                	  	((Q_temp.block_mat)[i*(Q_temp.n) + j]).n   = (a->n);
					}
                }
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
			Q.mat = malloc(Q.m * Q.n * sizeof(double));
    	    memcpy(Q.mat, mat_c.mat, Q.m * Q.n * sizeof(double));
			free(mat_c.mat);
	        rows_num = rows_num / 2;

			free_block_matrix(&Q_temp);
		}
    } // end of iter
	
	if (proc == 0){
        matrix q_temp;
        q_temp.m = loc_R.m;
        q_temp.n = loc_R.n;
        q_temp.mat = calloc( (q_temp.m) * (q_temp.n),  sizeof(double));
        
        qr_decomp(
            loc_R.m,
            loc_R.n,
            loc_R.mat,
            q_temp.mat,
            r->mat
        );
 
        mat_mul(&Q, &q_temp, q);
        free(q_temp.mat);
        free(Q.mat);
	}
}
