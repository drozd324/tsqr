//row major ideology
// electric fence for debugging!

#include "mat_tools.h"
#include "tsqr.h"

int main(){
    int m = 4;
    int n = 3;

    matrix a;
    a.mat = malloc(m * n * sizeof(double);
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
    
 

    return 0; 
}

