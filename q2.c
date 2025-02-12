//row major ideology
// electric fence for debugging!

#include "mat_tools.h"
#include "tsqr.h"

int main(){
	int m = 4;
	int n = 3;

	double* a = malloc(m * n * sizeof(double));
	srand(1234);
	for (int i=0; i<(m*n); i++){
		a[i] = (double)rand() / (double)RAND_MAX;
	}
	
	double* r = calloc(n * n, sizeof(double));

    return 0; 
}

