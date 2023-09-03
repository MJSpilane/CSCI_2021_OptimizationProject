// optimized versions of matrix diagonal summing
#include "matvec.h"

int matsquare_VER1(matrix_t *mat, matrix_t *matsq) {
    matrix_t matrix = *mat;
    matrix_t matsquared = *matsq;

    for(int i = 0; i < matrix.rows; i++){
        int j;
        for(j = 0; j < matrix.cols-4; j+=4){
            int k;
            //j + 0
            MSET(matsquared, i, j, 0);
            for(k=0; k<matrix.rows-4; k+=4){
                int mik = ((matrix).data[((i)*((matrix).cols)) + (k)]);
                int mkj = ((matrix).data[((k)*((matrix).cols)) + (j)]);
                
                int mik1 = ((matrix).data[((i)*((matrix).cols)) + (k+1)]);
                int mkj1 = ((matrix).data[((k+1)*((matrix).cols)) + (j)]);
                
                int mik2 = ((matrix).data[((i)*((matrix).cols)) + (k+2)]);
                int mkj2 = ((matrix).data[((k+2)*((matrix).cols)) + (j)]);
                
                int mik3 = ((matrix).data[((i)*((matrix).cols)) + (k+3)]);
                int mkj3 = ((matrix).data[((k+3)*((matrix).cols)) + (j)]);

                int val = mik*mkj;
                int val1 = mik1*mkj1;
                int val2 = mik2*mkj2;
                int val3 = mik3 *mkj3;
                
                int cur = MGET(matsquared, i, j);
                MSET(matsquared, i, j, (cur+val+val1+val2+val3));
            }
            for(; k<matrix.rows; k++){
                int mik = MGET(matrix, i, k);
                int mkj = MGET(matrix, k, j);
                int cur = MGET(matsquared, i, j);
                ((matsquared).data[((i)*((matsquared).cols)) + (j)] = (cur+(mik*mkj)));
            }
            //j+1
            ((matsquared).data[((i)*((matsquared).cols)) + (j+1)] = (0));
            for(k=0; k<matrix.rows-4; k+=4){
                int mik = ((matrix).data[((i)*((matrix).cols)) + (k)]);
                int mkj = ((matrix).data[((k)*((matrix).cols)) + (j+1)]);
                int val = mik*mkj;
                int mik1 = ((matrix).data[((i)*((matrix).cols)) + (k+1)]);
                int mkj1 = ((matrix).data[((k+1)*((matrix).cols)) + (j+1)]);
                int val1 = mik1*mkj1;
                int mik2 = ((matrix).data[((i)*((matrix).cols)) + (k+2)]);
                int mkj2 = ((matrix).data[((k+2)*((matrix).cols)) + (j+1)]);
                int val2 = mik2*mkj2;
                int mik3 = ((matrix).data[((i)*((matrix).cols)) + (k+3)]);
                int mkj3 = ((matrix).data[((k+3)*((matrix).cols)) + (j+1)]);
                int val3 = mik3 *mkj3;
                int cur = MGET(matsquared, i, j+1);
                ((matsquared).data[((i)*((matsquared).cols)) + (j+1)] = (cur + val + val1+val2+val3));

            }
            for(; k<matrix.rows; k++){
                int mik = MGET(matrix, i, k);
                int mkj = MGET(matrix, k, j+1);
                int cur = MGET(matsquared, i, j+1);
                ((matsquared).data[((i)*((matsquared).cols)) + (j+1)] = (cur + (mik*mkj)));

            }

            //j + 2
            ((matsquared).data[((i)*((matsquared).cols)) + (j+2)] = (0));
            for(k=0; k<matrix.rows-4; k+=4){
                int mik = ((matrix).data[((i)*((matrix).cols)) + (k)]);
                int mkj = ((matrix).data[((k)*((matrix).cols)) + (j+2)]);
                int val = mik*mkj;
                int mik1 = ((matrix).data[((i)*((matrix).cols)) + (k+1)]);
                int mkj1 = ((matrix).data[((k+1)*((matrix).cols)) + (j+2)]);
                int val1 = mik1*mkj1;
                int mik2 = ((matrix).data[((i)*((matrix).cols)) + (k+2)]);
                int mkj2 = ((matrix).data[((k+2)*((matrix).cols)) + (j+2)]);
                int val2 = mik2*mkj2;
                int mik3 = ((matrix).data[((i)*((matrix).cols)) + (k+3)]);
                int mkj3 = ((matrix).data[((k+3)*((matrix).cols)) + (j+2)]);
                int val3 = mik3 *mkj3;
                int cur = MGET(matsquared, i, j+2);
                ((matsquared).data[((i)*((matsquared).cols)) + (j+2)] = (cur + val + val1 + val2 + val3));
            }
            for(; k<matrix.rows; k++){
                int mik = MGET(matrix, i, k);
                int mkj = MGET(matrix, k, j+2);
                int cur = MGET(matsquared, i, j+2);
                ((matsquared).data[((i)*((matsquared).cols)) + (j+2)] = (cur + (mik*mkj)));

            }

            //j + 3 
            ((matsquared).data[((i)*((matsquared).cols)) + (j+3)] = (0));
            for(k=0; k<matrix.rows-4; k+=4){
                int mik = ((matrix).data[((i)*((matrix).cols)) + (k)]);
                int mkj = ((matrix).data[((k)*((matrix).cols)) + (j+3)]);
                int val = mik*mkj;
                int mik1 = ((matrix).data[((i)*((matrix).cols)) + (k+1)]);
                int mkj1 = ((matrix).data[((k+1)*((matrix).cols)) + (j+3)]);
                int val1 = mik1*mkj1;
                int mik2 = ((matrix).data[((i)*((matrix).cols)) + (k+2)]);
                int mkj2 = ((matrix).data[((k+2)*((matrix).cols)) + (j+3)]);
                int val2 = mik2*mkj2;
                int mik3 = ((matrix).data[((i)*((matrix).cols)) + (k+3)]);
                int mkj3 = ((matrix).data[((k+3)*((matrix).cols)) + (j+3)]);
                int val3 = mik3 *mkj3;
                int cur = MGET(matsquared, i, j+3);
                ((matsquared).data[((i)*((matsquared).cols)) + (j+3)] = (cur + val + val1 + val2 + val3));

            }
            for(; k<matrix.rows; k++){
                int mik = MGET(matrix, i, k);
                int mkj = MGET(matrix, k, j+3);
                int cur = MGET(matsquared, i, j+3);
                ((matsquared).data[((i)*((matsquared).cols)) + (j+3)] = (cur + (mik*mkj)));
            }
        }
        for(;j < matrix.cols;j++){
            int k;
            MSET(matsquared, i, j, 0);
            ((matsquared).data[((i)*((matsquared).cols)) + (j)] = (0));
            for(k=0; k<matrix.rows-4; k+=4){
                int mik = ((matrix).data[((i)*((matrix).cols)) + (k)]);
                int mkj = ((matrix).data[((k)*((matrix).cols)) + (j)]);
                int val = mik*mkj;
                int mik1 = ((matrix).data[((i)*((matrix).cols)) + (k+1)]);
                int mkj1 = ((matrix).data[((k+1)*((matrix).cols)) + (j)]);
                int val1 = mik1*mkj1;
                int mik2 = ((matrix).data[((i)*((matrix).cols)) + (k+2)]);
                int mkj2 = ((matrix).data[((k+2)*((matrix).cols)) + (j)]);
                int val2 = mik2*mkj2;
                int mik3 = ((matrix).data[((i)*((matrix).cols)) + (k+3)]);
                int mkj3 = ((matrix).data[((k+3)*((matrix).cols)) + (j)]);
                int val3 = mik3 *mkj3;
                int cur = MGET(matsquared, i, j);
                ((matsquared).data[((i)*((matsquared).cols)) + (j)] = (cur + val + val1 + val2 + val3));
            }
            for(; k<matrix.rows; k++){
                int mik = MGET(matrix, i, k);
                int mkj = MGET(matrix, k, j);
                int cur = MGET(matsquared, i, j);
                ((matsquared).data[((i)*((matsquared).cols)) + (j)] = (cur+ (mik*mkj)));
            }
        }
    }
    *matsq = matsquared;
    return 0;
}

int matsquare_VER2(matrix_t *mat, matrix_t *matsq) {
  // OPTIONALLY, OTHER VERSIONS
  matrix_t matrix = *mat;
  matrix_t matsquared = *matsq;
  matrix_t transposed;
  matrix_init(&transposed, matrix.rows, matrix.cols);
  for (int i = 0; i < matrix.rows; i++){
    for(int j = 0; j < matrix.cols; j++){

        MSET(transposed, i,j, MGET(matrix, j,i));

    }
  }

  for(int i = 0; i<matrix.rows; i++){
    for(int j = 0; j <matrix.cols; j++){
        MSET(matsquared, i, j, 0);
        int k;
        for(k = 0; k < matrix.rows-4; k+=4){
            //k+0
            int mik = ((matrix).data[((i)*((matrix).cols)) + (k)]);
            int mkj = ((transposed).data[((j)*((matrix).cols)) + (k)]);
            int val = mik*mkj;
            //k+1
            int mik1 = ((matrix).data[((i)*((matrix).cols)) + (k+1)]);
            int mkj1 = ((transposed).data[((j)*((matrix).cols)) + (k+1)]);
            int val1 = mik1*mkj1;

            //k+2
            int mik2 = ((matrix).data[((i)*((matrix).cols)) + (k+2)]);
            int mkj2 = ((transposed).data[((j)*((matrix).cols)) + (k+2)]);
            int val2 = mik2*mkj2;

            //k+3
            int mik3 = ((matrix).data[((i)*((matrix).cols)) + (k+3)]);
            int mkj3 = ((transposed).data[((j)*((matrix).cols)) + (k+3)]);
            int val3 = mik3*mkj3;

            int cur = ((matsquared).data[((i)*((matrix).cols)) + (j)]);
            int new = cur + val + val1 + val2 + val3;
            MSET(matsquared, i, j, new);
        }
        for(; k<matrix.rows; k++){
            int mik = ((matrix).data[((i)*((matrix).cols)) + (k)]);
            int mkj = ((transposed).data[((j)*((matrix).cols)) + (k)]);
            int cur = MGET(matsquared, i, j);
            int new = cur + mik* mkj;
            MSET(matsquared,i ,j, new);
        }
    }
  }
  matrix_free_data(&transposed);
  *matsq = matsquared;
  return 0;
}


int matsquare_OPTM(matrix_t *mat, matrix_t *matsq){
  if(mat->rows != mat->cols   ||      // must be a square matrix to square it
     mat->rows != matsq->rows ||
     mat->cols != matsq->cols)
  {
    printf("matsquare_OPTM: dimension mismatch\n");
    return 1;
  }

  // Call to some version of optimized code
  return matsquare_VER2(mat, matsq);
}



/////////////////////////////////////////////////////////////////////////////////
// ADDITIONAL INFORMATION
//
// (A) VERSION: If you implemented several versions, indicate which
// version you timed
// 
// ####################### YOUR ANSWER HERE #########################
// VER2
// ##################################################################
// 
//
// (B) TIMING ON loginNN.cselabs.umn.edu:
// Paste a copy of the results of running matsquare_benchmark on the
// above machines in the space below which shows how your performance
// optimizations improved on the baseline codes.
// 
// ####################### YOUR ANSWER HERE #########################
//==== Matrix Square Benchmark Version 2 ====
//  SIZE       BASE       OPTM  SPDUP   LOG2  SCALE POINTS 
//   256 3.8801e-01 4.2758e-02   9.07   3.18   0.94   2.98 
//   273 4.0617e-01 5.1489e-02   7.89   2.98   1.00   2.98 
//   512 3.9508e+00 3.4608e-01  11.42   3.51   1.88   6.59 
//   801 1.4735e+01 1.3711e+00  10.75   3.43   2.93  10.05 
//  1024 5.3791e+01 2.8468e+00  18.89   4.24   3.75  15.90 
//RAW POINTS: 38.51
//TOTAL POINTS: 30 / 30
// ##################################################################
// 
// (C) OPTIMIZATIONS:
// Describe in some detail the optimizations you used to speed the code
// up.  THE CODE SHOULD CONTAIN SOME COMMENTS already to describe these
// but in the section below, describe in English the techniques you used
// to make the code run faster.  Format your descriptions into discrete
// chunks such as.
// 
// Optimization 1: Blah bla blah... This should make run faster because
// yakkety yakeety yak.
// 
// Optimization 2: Blah bla blah... This should make run faster because
// yakkety yakeety yak.
// ...
// Optimization N: Blah bla blah... This should make run faster because
// yakkety yakeety yak.
// 
// Full credit solutions will describe 2-3 optimizations and describe
// WHY these improved performance in at least a couple sentences.
// 
// ####################### YOUR ANSWER HERE #########################
//Optimization 1: Local Copies of Variables. Creating local copies of the matrices improves 
// performance because it will create copies of this data and store it in memory that is accessed
// much more quickly than the main memory wher the matrices are located.
//
//Optimization 2: Replace function calls. Function calls are slow because on the assembly level
//they require jumping and restrict register usage. Replacing these with Macros that are inserted
//directly into the code at the assembly level improves performance, however, simply putting the macro
//code in where the macros would go improves performance even more at the cost of readability.
//
// Optimization 3: Transposing the matrix. Even though it adds two more for loops,
// these for loops only iterate once and it allows for the transposed matrix to be
// iterated row-wise and give equivalent values as iterating the original matrix
// column-wise, because C is a row major language, this will greatly increase performance
// due to many less cache misses. 
//
//Optimization 4: Loop Unrolling. Loop unrolling was utilized on the inner most loop
//(the k loop). This improves performance because it allows super scalar CPUs to utilize multiple 
// arithmetic  logic units to do multiple computations at a time, rather than doing one value of k
//each loop, which forces the CPU to do computations much more linearly, and thus, slowly.
// ##################################################################
