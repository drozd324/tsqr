//row major ideology
// electric fence for debugging!

//#include "mat_tools.h"
#include "tsqr.h"
#include <cblas.h>
#include <stdlib.h>

int main(){
    int m = 12;
    int n = 4;

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
    r.m = m;
    r.n = n;

    tsqr(&a, &q, &r); 
	
    matrix qr;
    qr.mat = calloc(m * n, sizeof(double));
    qr.m = m;
    qr.n = n;

    mat_mul(&q, &r, &qr);

    printf("printing input matrix a");
    print_matrix(a.m, a.n, a.mat, a.n);

    printf("printing input matrix q");
    print_matrix(q.m, q.n, q.mat, q.n);

    printf("printing input matrix r");
    print_matrix(r.m, r.n, r.mat, r.n);

    printf("printing input matrix qr");
    print_matrix(qr.m, qr.n, qr.mat, qr.n);

    return 0; 
}

