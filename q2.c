#include <stdio.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>

//row major ideology
// electric fence for debugging!

typedef struct _Matrix{
	double* mat;
	int m;
	int n;
}matrix;

typedef struct _block_Matrix{
	matrix* block_mat;
	int m;
	int n;
}block_matrix;

void print_matrix(int m, int n, double* a, int lda);
void qr_decomp(int m, int n, double* a, double* q, double* r);
void decomp1d(int n, int N, int block_index, int s, int e);

void decomp_matrix(matrix* a, block_matrix* block_a);
void free_block_matrix(block_matrix* block_a);
void comp_matrix(block_matrix* block_a, matrix* a);
void tsqr(matrix* a, matrix* q, matrix* r);

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

/**
 * @brief Prints a matrix from a single contigious array
 *
 * @param m Number of rows of matrix
 * @param n Number of columns of matrix
 * @param a Pointer to array of doubles representing matrix
 * @param lda Leading dimension of the matrix 
*/
void print_matrix(int m, int n, double* a, int lda) {
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++){
            printf("%2.2f ", a[i*lda + j]);
		}
        printf("\n");
    }
}

/**
 * @brief A function which combines two of LAPACKE's qr factorisation
 *        function to obain 
 *
 * @param m Number of rows of matrix to decompose
 * @param n Number of columns of matrix to decompose
 * @param a Pointer to array of doubles representing matrix to decompose
 * @param q Pointer to array of doubles representing matrix to put Q
 *          from QR = A into.
 * @param r Pointer to array of doubles representing matrix to put R
 *          from QR = A into.
*/
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


/**
 * @brief Takes in the lenght of an array n, the desired amount of blocks 
 *        to split it into N, the index of the desired block block_index,
 *        and gives the user start and endpoints for this block int the 
 *        original array s, e.
 * 
 * @param n Number of entries in array 
 * @param N Number of desired blocks 
 * @param block_index What would be the index if the block 
 * @param s Start of block indexed by block_index in original array
 * @param e End of block indexed by block_index in original array
*/
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

/**
 * @brief Decomposes an input matrix struct into a block_matrix struct
 *        with dimensions of block matrix given in the block_matrix 
 *        struct.
 * 
 * @param a Pointer to an empty matrix struct
 * @param block_a Pointer to a block_matrix struct with assigned shape, (n and m).
*/
void decomp_matrix(matrix* a, block_matrix* block_a){
	int s1;
	int e1;
	int s2;
	int e2;
    
    // for each desired block
	for (int i=0; i < (block_a->m); i++){
        decomp1d(a->m, block_a->m, i, s1, e1);
		for (int j=0; j < (block_a->n); j++){
            decomp1d(a->n, block_a->n , j, s2, e2);
           
            // make a block
			matrix block;
			block.m = e1 - s1;
			block.n = e2 - s2;
			block.mat = malloc((block.m) * (block.n) * sizeof(double)); 
			// gotta free this memory yourself after 
			// you're finished with the block matrix
            
            // iterate over each entry of this block and plop in corresponding 
            // values from the matrix struct a
			for (int k=0; k < block.m; k++){
				for (int l=0; l < block.n; l++){
                    (block.mat)[k*(block.n) + l] = (a->mat)[(s1+k)*(a->n) + s2+l];
				}
			}	
            
            // assign address of this new martix struct "block" to an 
            // entry in the new block matrix
			(block_a->block_mat)[i*(block_a->n) + j] = &block
		
		}
	}
}


/**
 * @brief Frees allocated memory by the decomp_matrix function. 
 * 
 * @param block_a Pointer to a block_matrix struct with assigned shape, (n and m).
*/
void free_block_matrix(block_matrix* block_a){
    // iterate over each entry of block matrix
	for (int i=0; i<block_a->m; i++){
		for (int j=0; j<block_a->n; j++){
            
            // picks out address of array within matrix struct, within block matrix to free
			free(&( ( (block_a->block_mat)[i*(block_a->n) + j] )->mat ) );
}


//void comp_matrix(int m, int n, double* a, int M, int N, double* a){
/**
 * @brief Reconstructs (composes) a matrix struct back from a block_block matrix struct
 *        created by decomp_matrix.
 * 
 * @param block_a Pointer to a block_matrix struct with assigned shape, (n and m).
 * @param a Pointer to an empty matrix struct
*/
void comp_matrix(block_matrix* block_a, matrix* a){
    int row_step = 0;
    int col_step = 0;
    int sub_rows = 0;
    int sub_cols = 0;
    
     
	for (int i=0; i<block_a->m; i++){
		for (int j=0; j<block_a->n; j++){
            
			sub_rows = ((block_a->block_mat)[i*(block_a->n) + j])->m;
			sub_cols = ((block_a->block_mat)[i*(block_a->n) + j])->n;

        	for (int k=0; k<sub_rows; k++){
        		for (int l=0; l<sub_rows; l++){
                    a.mat[(row_step + k)*n + (col_step + l)] = block_a[k*sub_cols + l];

            // picks out address of array within matrix struct, within block matrix to free
			free(&( ( (block_a->block_mat)[i*(block_a->n) + j] )->mat ) );

            col_step += sub_cols;
        row_step += sub_rows;
        col_step = 0;
	}
}

// calloc q.matrix
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
		q.mat[i*n + i] = 1;
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
		Q_temp->block_mat = malloc(rows_num * rows_num * sizeof(matrix*));
        
        // init R		
		block_matrix* R_temp;
        R_temp->m = rows_num;
        R_temp->n = 1;
		R_temp->block_mat = malloc(rows_num * sizeof(matrix*));
        
        // allocate memory for Q_temp
        for (int i=0; i<rows_num; i++){
            for (int j=0; j<rows_num; j++){
                // (array of matrix structs)[of matrix struct]
                // (block matrix)[index of matrix]->(array of double representing matrix)  =  calloc size m_i*n array of doubles representing matrix
                ((Q_temp->block_mat)[i*(Q_temp->n) + j])->mat = calloc(((R->block_mat)->m) * (a->n), sizeof(double*));
                ((Q_temp->block_mat)[i*(Q_temp->n) + j])->m   = (R->block_mat)->m;
                ((Q_temp->block_mat)[i*(Q_temp->n) + j])->n   = a->n;
            }
        } 
        
        // allocate memory for R_temp
        for (int i=0; i<rows_num; i++){
            ((R_temp->block_mat)[i])->mat = calloc((a->n) * (a->n), sizeof(double*));
            ((R_temp->block_mat)[i])->m   = (R->block_mat)->m;
            ((R_temp->block_mat)[i])->n   = a->n;
        }
		
		for (int i=0; i<rows_num; i++){
            // index the block matrix and pick out the array of doubles representing matrix 
            //        m  n  (mat to decomp)
			qr_decomp(m, n, ((R->block_mat)[i])->mat, ((Q_temp->block_mat)[i*(Q_temp->n) + i])->mat, ((R_temp->block_mat)[i])->mat;
		    //                                      (          pick out  (i,i) entry in Q_temp   ) (pick out (i,0) entry in block mat)
		    //                                      (          which would be a matrix           ) ( which is again a matrix         )
        } 
        
        // matrix needed for computation
        matrix* Q_temp_matrix;
        Q_temp_matrix->m = (rows_num == max_block_rows) ? m : rows_num * 2 * (a->n);
        Q_temp_matrix->n = rows_num * (a->n);
        comp_matrix(Q_temp, Q_temp_matrix); // makes the block_matrix struct Q_temp into a matrix struct Q_temp_matrix
	 
        // stuff for matrix mul to follow
        int k = fmin(a->m, Q_temp_matrix->n);
        int alpha = 1;
        int beta = 0;
        int lda = a->n;
 
 	    double* c = calloc(m * n, sizeof(double));
 	    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
 	   				a->m, Q_temp_matrix->n, k, alpha, q->mat, lda,
 	   				Q_temp_comp->mat, lda, beta,
 	   				c, lda);
                    
	    memcpy(q->mat, c, n*m * sizeof(double));
    
        rows_num = rows_num / 2;
        free_block_mat(R);
        R->m = rows_num;
        for (int i=0; i<rows_num; i++){
            block_mat temp_block;
            temp_block.
            
    
            (R->block_mat)[]

        }
        
        
                     
                
        
        




        
	}
    
}
