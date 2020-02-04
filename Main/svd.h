#pragma once
#include <memory>

struct SVD {
  const double *U;
  const double *s;
  const double *VT;
  const int k;
  const int m;
  const int n;

      ~SVD() {
        free(const_cast<double *>(U));
        free(const_cast<double *>(s));
        free(const_cast<double *>(VT));
    }

};

void free_svd_result(SVD* svd);
void transpose( const double * __restrict__  input, double * __restrict__ output, const int N, const int M);
void matrix_matrixTranspose_Multi(const double * __restrict__  input_1 , double * __restrict__ output, const int n, const int m,  const int SIZE_Out);

SVD* compute_svd(double *A, int m, int n,
                        double threshold = -1.0,
                        double accuracy = 1e-6);