//row major ideology
// electric fence for debugging!

//#include "mat_tools.h"
#include "tsqr.h"
#include <cblas.h>
#include <lapacke.h>
#include <stdlib.h>

int main(){
    int m = 20;
    int n = 5;

    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("~~~~~~~~~~~~~~~~~~~~             q2.c                ~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

    for (int func_type=0; func_type<2; func_type++){

        matrix a;
        a.mat = malloc(m * n * sizeof(double));
        a.m = m;
        a.n = n;
        srand(1234);
        for (int i=0; i<(m*n); i++){
            (a.mat)[i] = (double)rand() / (double)RAND_MAX;
        }
       
        matrix q; 
        q.mat = calloc(m * n, sizeof(double));
        q.m = m;
        q.n = n;
 
        matrix r; 
        r.mat = calloc(n * n, sizeof(double));
        r.m = n;
        r.n = n;
	
	if (func_type==0){
	    printf("-------------------- Using TSQR -----------------------\n");
            tsqr(&a, &q, &r); 
	} else {
	    printf("------------------ Using serial QR --------------------\n");
            qr_decomp(a.m, a.n, a.mat, q.mat, r.mat);
        } 
            
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
    return 0; 
}

