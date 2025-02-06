#include <stdio.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>

//row major ideology
//change struct to struct pointers in functions pass

struct _Matrix{
	double* mat;
	int m;
	int n;
}matrix;

struct _block_Matrix{
	matrix* block_mat;
	int m;
	int n;
}block_matrix;

void print_matrix(int m, int n, double* a, int lda); 
void qr_decomp(int m, int n, double* a, double* r);
void decomp1d(n, N, block_index);

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


void print_matrix(int m, int n, double* a, int lda) {
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++){
            printf("%2.2f ", a[i*lda + j]);
		}
        printf("\n");
    }
}

void qr_decomp(int m, int n, double* a, double* q, double* r){
    int     lda = n;
	double* tau = malloc(fmin(m, n) * sizeof(double));
	
	memcpy(q, a, n*m * sizeof(double));
	
	// Double precision GEneral QR Factorisation
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, m, n, q, lda, tau); 
	
	// write to R to r
    for (int i=0; i<n; i++) {
        for (int j=i; j<n; j++){
            r[i*n + j] = q[i*n + j];
		}
	}

	int k = fmin(m, n);
	// Double precision ORthogonal Generate? QR Factorisation
   	LAPACKE_dorgqr(LAPACK_ROW_MAJOR, m, n, k, q, lda, tau); 
}


void decomp1d(int n, int N, int block_index, int s, int e){
    int remainder = n % N;
    int base = n / N;

    if (block_index < remainder){
        s =  block_index * (base+1);
        e = s + base+1;
    } else {
        s = (remainder * (base+1)) + ((fmax(block_index-remainder, 0)) * base);
        e = s + base;
    }
}

//void decomp_matrix(matrix* a, int M, int N, block_matrix* block_a){
//block_a->m = M;
//block_a->n = N;	
//block_a->block_matrix = malloc(block->m * block->n * sizeof(matrix*));
void decomp_matrix(matrix* a, block_matrix* block_a){
	int s1;
	int e1;
	int s2;
	int e2;

	for (int i=0; i < (block_a->m); i++){
        decomp1d(a.m, block_a->m, i, s1, e1);
		for (int j=0; j < (block_a->n); j++){
            decomp1d(a.n, block_a->n , j, s2, e2);
           
			matrix block;
			block.m = e1 - s1;
			block.n = e2 - s2;
			
			block.mat = malloc((block.m) * (block.n) * sizeof(double)); 
			// gotta free this memory yourself after 
			// you're finished with the block matrix
			for (int k=0; k < block.m; k++){
				for (int l=0; l < block.n; l++){
                    (block.mat)[k*(block.n) + l] = (a->mat)[(s1+k)*(a->n) + s2+l];
				}
			}	

			(block_a->block_mat)[i*n + j] = &block
		
		}
	}
}



//void comp_matrix(int m, int n, double* a, int M, int N, double* a){
void comp_matrix(block_matrix* block_a, matrix* a){
    int row_step = 0;
    int col_step = 0;
    int sub_rows = 0;
    int sub_cols = 0;
    
	for (int i=0; i<block_a->m; i++){
		for (int j=0; j<block_a->n; j++){
            
			sub_rows = (block_a.block_mat[i*N + j])->m;
			sub_cols = (block_a.block_mat[i*N + j])->n;

        	for (int k=0; k<sub_rows; k++){
        		for (int l=0; l<sub_rows; l++){
                    a.mat[(row_step + k)*n + (col_step + l)] = block_a[k*sub_cols + l];

			// maybe
			free(&(block_a.block_mat[i*N + j]));

            col_step += sub_cols;
        row_step += sub_rows;
        col_step = 0;
	}
}

// calloc q.matrix
void tsqr(matrix* a, matrix* q, matrix* r){
	for (int i=0; i<m; i++){
		q.mat[i*n + i] = 1;
	}

	int rows_num = 4;
	block_matrix* A;
    A->m = rows_num;
    A->n = 1;
	decomp_matrix(a, A);
	
	for (int iter=0; iter<3; iter++){
		
		block_matrix* Q_temp;
        Q_temp->m = rows_num;
        Q_temp->n = rows_num;
		Q_temp->block_mat = malloc(rows_num * rows_num * sizeof(matrix*));
		
		block_matrix* R_temp;
        R_temp->m = rows_num;
        R_temp->n = 1;
		R_temp->block_mat = malloc(rows_num * sizeof(matrix*));
        
        // allocate memory for Q_temp
        for (int i=0; i<rows_num; i++){
            for (int j=0; j<rows_num; j++){
                // (array of matrix structs)[of matrix struct]
                // (block matrix)[index of matrix]->(array of double representing matrix)  =  calloc size m_i*n array of doubles representing matrix
                (Q_temp->block_mat)[i*(Q_temp->n) + j]->mat = calloc((A->block_mat->m) * (a->n), sizeof(double*));
                (Q_temp->block_mat)[i*(Q_temp->n) + j]->m   = A->block_mat->m;
                (Q_temp->block_mat)[i*(Q_temp->n) + j]->n   = a->n;
            }
        } 
        
        // allocate memory for R_temp
        for (int i=0; i<rows_num; i++){
            (R_temp->block_mat)[i]->mat = calloc((a->n) * (a->n), sizeof(double*));
            (R_temp->block_mat)[i]->m   = A->block_mat->m;
            (R_temp->block_mat)[i]->n   = a->n;

        }
		
		for (int i=0; i<rows_num; i++){
			qr_decomp(m, n, R->mat,
                // index the block matrix and pick out the array of doubles representing matrix 
                (Q_temp->block_mat)[i*(Q_temp->n) + i]->mat, 
                (R_temp->block_mat)[i                ]->mat);
		} 
	
        for (int i=0; i<rows_num; i++){
            for (int j=0; j<rows_num; j++){
                if (i != j){
                    (Q_temp->block_mat)[i*(Q_temp->n) + j] 
                    = calloc((Q_temp->m) * (Q_temp->n), sizeof(double)); 
                }
            }
        }	

	}
	

}
