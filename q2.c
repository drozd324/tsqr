//row major ideology
// electric fence for debugging!

//#include "mat_tools.h"
#include "tsqr.h"
#include <cblas.h>
#include <lapacke.h>
#include <stdlib.h>
#include <mpi.h>

int main(){
	int rank;
	MPI_Init(NULL, NULL);		
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);	

    matrix a;
    matrix q;
    matrix r;
    int m = 20;
    int n = 5;
	if (rank == 0){
 
        printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        printf("~~~~~~~~~~~~~~~~~~~~             q2.c                ~~~~~~~~~~~~~~~~~~~~~~~~\n");
        printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
 
 
        a.mat = malloc(m * n * sizeof(double));
        a.m = m;
        a.n = n;
        srand(1234);
        for (int i=0; i<(m*n); i++){
            (a.mat)[i] = (double)rand() / (double)RAND_MAX;
        }
        
        q.mat = calloc(m * n, sizeof(double));
        q.m = m;
        q.n = n;
  
        r.mat = calloc(n * n, sizeof(double));
        r.m = n;
        r.n = n;
	
		printf("-------------------- Using TSQR -----------------------\n");
 	}
 
    tsqr(&a, &q, &r); 
 	   	
    
	if (rank == 0){        
        matrix qr;
        qr.mat = calloc(m * n, sizeof(double));
        qr.m = m;
        qr.n = n;
  
        mat_mul(&q, &r, &qr);
  
        printf("printing input matrix a\n");
        print_matrix(a.m, a.n, a.mat, a.n);
        printf("\n");
  
        printf("printing output matrix q\n");
        print_matrix(q.m, q.n, q.mat, q.n);
        printf("\n");
  
        printf("printing output matrix r\n");
        print_matrix(r.m, r.n, r.mat, r.n);
        printf("\n");
  
        printf("printing output matrix qr\n");
        print_matrix(qr.m, qr.n, qr.mat, qr.n);
        printf("\n");
  
		printf("------------------ Using serial QR --------------------\n");
		qr_decomp(a.m, a.n, a.mat, q.mat, r.mat);
        mat_mul(&q, &r, &qr);
            
        printf("printing input matrix a\n");
        print_matrix(a.m, a.n, a.mat, a.n);
        printf("\n");
  
        printf("printing output matrix q\n");
        print_matrix(q.m, q.n, q.mat, q.n);
        printf("\n");
  
        printf("printing output matrix r\n");
        print_matrix(r.m, r.n, r.mat, r.n);
        printf("\n");
  
        printf("printing output matrix qr\n");
        print_matrix(qr.m, qr.n, qr.mat, qr.n);




        printf("\n");
        free(a.mat);
        free(q.mat);
        free(r.mat);
        free(qr.mat);


        /*
        double avg = 0;
        for (int i=0; i<m*n; i++){
            avg += (a.mat[i] - qr.mat[i]) * (a.mat[i] - qr.mat[i]);
        }
        avg = avg / (m*n);
  
        printf("avg = %lf\n", avg);
        */
    }

	MPI_Finalize();
    return 0; 
}

