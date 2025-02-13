/**
 * @file mat_tools.h 
 * @brief Header for mat_tools.c 
 *
 * @author Patryk Drozd
 */


/**
 * @brief Struct for matrices
*/
typedef struct _matrix{
	double* mat;  // single coigious array of doubles for entries of matrix
	int m;  // number of rows of matrix
	int n;  // number of columns of matrix
}matrix;

/**
 * @brief Struct for block matrices
*/
typedef struct _block_matrix{
	matrix* block_mat;  // single contigious array for pointers to matrix structs
	int m;  // number of rows of blocks of block matrix
	int n;  // number of columns of block of block matrix
}block_matrix;


void print_matrix(int m, int n, double* a, int lda);
void qr_decomp(int m, int n, double* a, double* q, double* r);
void mat_mul(matrix* a, matrix* b, matrix* c);
void decomp1d(int n, int N, int block_index, int* s, int* e);

void decomp_matrix(matrix* a, block_matrix* block_a);
void free_block_matrix(block_matrix* block_a);
void comp_matrix(block_matrix* block_a, matrix* a);

