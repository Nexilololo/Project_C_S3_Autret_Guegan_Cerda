#include "matrix.h"

t_matrix createZeroMatrix(int n){
    t_matrix M;
    M.rows = n;
    M.cols = n;
    M.data = (double**)calloc(n, sizeof(double*));
    for(int i = 0; i < n; i++){
        M.data[i] = (double*)calloc(n, sizeof(double));
    }
    return M;
}

void freeMatrix(t_matrix *M) {
    if (M->data != NULL) {
        for (int i = 0; i < M->rows; i++) {
            free(M->data[i]);
        }
        free(M->data);
        M->data = NULL;
    }
}

t_matrix createTransitionMatrix(t_adjlist graph){
    int n = graph.list_number;
    t_matrix M = createZeroMatrix(n);
    for(int i = 0; i < n; i++){
        t_cell *curr = graph.list[i].head;
        while(curr != NULL){
            int j = curr->vertex - 1;
            M.data[i][j] = curr->weight;
            curr = curr->next;
        }
    }
    return M;
}

void copyMatrix(t_matrix *dest, t_matrix src){
    for(int i = 0; i < src.rows; i++){
        for(int j = 0; j < src.cols; j++){
            dest->data[i][j] = src.data[i][j];
        }
    }
}

t_matrix multiplyMatrix(t_matrix A, t_matrix B) {
    // 1. Validation: Columns of A must equal Rows of B
    if (A.cols != B.rows) {
        printf("ERROR: Dimension mismatch! A: %dx%d, B: %dx%d\n", A.rows, A.cols, B.rows, B.cols);
        exit(1);
    }

    // 2. Create Result Matrix of correct size (A.rows x B.cols)
    t_matrix result;
    result.rows = A.rows;
    result.cols = B.cols;

    // Manually allocate memory because createZeroMatrix might be square-only
    result.data = (double**)calloc(result.rows, sizeof(double*));
    for (int i = 0; i < result.rows; i++) {
        result.data[i] = (double*)calloc(result.cols, sizeof(double));
    }

    // 3. Perform Multiplication
    for (int i = 0; i < result.rows; i++) {       // Rows of A
        for (int j = 0; j < result.cols; j++) {   // Cols of B
            double sum = 0.0;
            for (int k = 0; k < A.cols; k++) {    // Common dimension
                sum += A.data[i][k] * B.data[k][j];
            }
            result.data[i][j] = sum;
        }
    }

    return result;
}

double diffMatrix(t_matrix A, t_matrix B){
    double result = 0.0;
    for(int i = 0; i < A.rows; i++){
        for(int j = 0; j < A.cols; j++){
            result += fabs(A.data[i][j] - B.data[i][j]);
        }
    }
    return result;
}

void printMatrix(t_matrix M){
    for(int i = 0; i < M.rows; i++){
        for(int j = 0; j < M.cols; j++){
            printf("%.4f ", M.data[i][j]);
        }
        printf("\n");
    }
}

t_matrix subMatrix(t_matrix matrix, t_partition part, int compo_index){
    t_class currclass = part.partition[compo_index];
    int n = currclass.vertex_number;
    t_matrix sub = createZeroMatrix(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int original_row_index = currclass.vertex[i].ID - 1;
            int original_col_index = currclass.vertex[j].ID - 1;

            sub.data[i][j] = matrix.data[original_row_index][original_col_index];
        }
    }
    return sub;
}