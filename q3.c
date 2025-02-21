//#include "mat_tools.h"
#include "tsqr.h"
#include <cblas.h>
#include <lapacke.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <sys/time.h>


double walltime(){
	struct timeval t;
	gettimeofday(&t, NULL);
	double wtime = (double) (t.tv_sec + t.tv_usec*1e-6);
	return wtime;
}	


int main(){
	int rank;
	MPI_Init(NULL, NULL);		
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);	

    matrix a;
    matrix q;
    matrix r;
	
    int m;
    int n;
    srand(1234);
	int num_matrices = 100;
	int min_n = 10;
	int max_size = 500;
	double t0;
	FILE *pf;

	if (rank == 0){
		pf = fopen("q3.txt", "w");
		
    	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~  q3.c  ~~~~~~~~~~~~~~~~~~~~~~~~\n");
    	printf("Generating %d matrices with min(n)=%d and largest dimension = %d \n", num_matrices, min_n, max_size);
	}

	for (int i=0; i<num_matrices; i++){
		if (rank == 0){
			fprintf(stderr, "matrix %d/%d \r", i, num_matrices);
			
    		n = min_n + (rand() % max_size);
	    	m = n*10 ; // makes sure the matrix will be "tall and skinny"

            a.mat = malloc(m * n * sizeof(double));
            a.m = m;
            a.n = n;
            for (int i=0; i<(m*n); i++){
                (a.mat)[i] = (double)rand() / (double)RAND_MAX;
            }
            
            q.mat = calloc(m * n, sizeof(double));
            q.m = m;
            q.n = n;
      
            r.mat = calloc(n * n, sizeof(double));
            r.m = n;
            r.n = n;
			
			t0 = walltime();
		
 	 	}
	
    	tsqr(&a, &q, &r); 
		 	   
		if (rank == 0){
    		fprintf(pf, "%d %d %lf \n", m, n, walltime() - t0);

       		free(a.mat);
        	free(q.mat);
        	free(r.mat);
		}	 
    }
	
	if (rank == 0){		
		fclose(pf);
		printf("\n");
	}


	MPI_Finalize();
    return 0; 
}

