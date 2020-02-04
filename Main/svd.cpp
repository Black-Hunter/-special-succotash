#include "svd.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cstring>
#include "../headers/aligned_allocator.h"
#include <memory>
#include <thread>

template <class T>
using aligned_vector = std::vector<T, alligned_allocator<T, 64>>;
using namespace std;

// in AT n = m , and m = n
void matrix_matrixTranspose_Multi(const double * __restrict__  input_1 , double * __restrict__ output, const int n, const int m,  const int SIZE_Out){

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            #pragma omp simd aligned(input_1, output: 64)
            for (int k = 0; k < m; ++k) {
                output[i*SIZE_Out +j] += input_1[i*m + k] * input_1[j*m + k];
            }
        }
    }
}

void transpose( const double * __restrict__  input, double * __restrict__ output, const int N, const int M) {

#pragma omp /*parallel for */simd aligned(input, output: 64)
    for(int n = 0; n<N*M; n++) {
        int i = n/N;
        int j = n%N;
        output[n] = input[M*j + i];
    }
}

inline void eigenvectors(
        double * __restrict__ const  inMatrix,
        double * __restrict__ Out_S,
        double * __restrict__ U,
        const int k,
        int & counter,
        const double threshold,
        const double accuracy) {


    //cout<<"Thread 2\n";
    double lambda = 0.0;
    counter = 0;
    aligned_vector<double> tmp(k);

    //#pragma omp simd
    for (int i = 0; i < k; ++i)
        tmp.at(i) = 1.0;
    aligned_vector <double>    A_tmp(k*k) ;//= (double *) malloc(k * k* sizeof(double));
    memcpy(&(*A_tmp.data()), inMatrix, k * k * sizeof(double));

    for (int i = 0; i < k; i++) {

        double fact[] = {-1, 0.0};
        double deffrance = INT64_MAX;

        while (deffrance > accuracy) {

            fact[0] = lambda;
            aligned_vector<double> out(k);
            double * _out= out.data();
            //#pragma omp parallel for
            for (int a = 0; a < k; a++) {
                _out[a] = 0;
                #pragma omp simd aligned(_out: 64)
                for (int j = 0; j < k; j++) {
                    _out[a] += (A_tmp[a * k + j] * tmp.at(j));
                }
            }

            memcpy(&(*tmp.begin()), &(*out.begin()), k * sizeof(double));
            lambda = 0.0;

            #pragma omp simd reduction(+ : lambda)
            for (int a = 0; a < k; a++)
                lambda = lambda + (tmp.at(a) * tmp.at(a));

            lambda = sqrt(lambda);
            fact[1] = lambda;

            if (lambda != 0.0)
                for (int a = 0; a < k; ++a)
                    tmp.at(a) = tmp.at(a) * (1.0 / lambda);
            else
                break;

            deffrance = abs(fact[0] - fact[1]);
        }

        double lambda_sing = sqrt(lambda);

        if (threshold && lambda_sing < threshold)
            break;
        else
            counter++;

        Out_S[i] = lambda_sing;

        for (int j = 0; j < k; ++j) {
            U[i * k + j] = tmp.at(j);
        }

        aligned_vector<double> tmp_(k * k);
        for (int j = 0; j < k; ++j) {
            for (int l = 0; l < k; ++l) {
                tmp_.at(j * k + l) = (tmp.at(j) * tmp.at(l));
            }
        }

        double * _tmp_ = tmp_.data();
        #pragma omp simd aligned(_tmp_: 64)
        for (int j = 0; j < k * k; ++j) {
            _tmp_[j] = _tmp_[j] * lambda;
        }

        for (int j = 0; j < k * k; ++j) {
            A_tmp[j] = abs(A_tmp[j] - tmp_.at(j));
        }
    }

}

SVD * power_methode_deflation(const double *  A , const int n, const int m , const double accuracy, double threshold = 0.0) {
/*
    double * UT = (double *) malloc(n*n*sizeof(double));
    double * S  = (double *) malloc(max(n,m)*sizeof(double));
    double * V  = (double *) malloc(m*m*sizeof(double));
    double * AT = (double *) malloc(m*n*sizeof(double));
    double * ATA= (double *) malloc(m*m*sizeof(double));;
    double * AAT= (double *) malloc(n*n*sizeof(double));
 */
    #ifdef _MSC_VER
    const int N = 64;
    double * UT = (double *) _aligned_malloc(n*n*sizeof(double));
    double * S  = (double *) _aligned_malloc(max(n,m)*sizeof(double));
    double * V  = (double *) _aligned_malloc(m*m*sizeof(double));
    double * AT = (double *) _aligned_malloc(m*n*sizeof(double));
    double * ATA= (double *) _aligned_malloc(m*m*sizeof(double));;
    double * AAT= (double *) _aligned_malloc(n*n*sizeof(double));
    //(double *)_aligned_malloc(n * sizeof(value_type), N);
    #else
    const int N = 64;
    void * UT = nullptr;//(double *) malloc(n*n*sizeof(double));
    void * S  = nullptr;//(double *) malloc(max(n,m)*sizeof(double));
    void * V  = nullptr;//(double *) malloc(m*m*sizeof(double));
    void * AT = nullptr;//(double *) malloc(m*n*sizeof(double));
    void * ATA= nullptr;//(double *) malloc(m*m*sizeof(double));;
    void * AAT= nullptr;//(double *) malloc(n*n*sizeof(double));

    bool error = posix_memalign(&UT , N, n*n      * sizeof(double))||
                 posix_memalign(&S  , N, max(n,m) * sizeof(double))||
                 posix_memalign(&V  , N, m*m      * sizeof(double))||
                 posix_memalign(&AT , N, n*m      * sizeof(double))||
                 posix_memalign(&ATA, N, m*m      * sizeof(double))||
                 posix_memalign(&AAT, N, n*n      * sizeof(double));
        if (error) {
            cout<<"Error, couldn't allocate enough memory";
            exit(1); // couldn't allocate enough memory
        }
    #endif

    memset(UT, 0, n*n*sizeof(double));
    memset(V, 0, m*m*sizeof(double));
    memset(S, 0, max(n,m)*sizeof(double));

    int counter1 = 1;
    int counter2 = 1;

    // Create AT
    transpose((double *)A,(double *)AT, n,m);


        thread  th_matrix_matrixTranspose_Multi(matrix_matrixTranspose_Multi,(double *)AT, (double *)ATA, m, n, m);
        thread  th_matrix_matrixTranspose_Multi_(matrix_matrixTranspose_Multi,(double *)A, (double *)AAT, n, m, n);

        th_matrix_matrixTranspose_Multi.join();
        th_matrix_matrixTranspose_Multi_.join();

        thread th_claculateEigenvectors_(eigenvectors, (double *)ATA, (double *)S, (double *)V, m, ref(counter1), threshold, accuracy);
        thread th_claculateEigenvectors(eigenvectors, (double *)AAT, ((double *)S), (double *)UT, n, ref(counter2), threshold, accuracy);

        th_claculateEigenvectors_.join();
        th_claculateEigenvectors.join();

        int counter = min(counter1,counter2);


    S  =(double *) realloc(S , counter* sizeof(double));
    UT =(double *) realloc(UT, counter2*n* sizeof(double));
    V  =(double *) realloc(V , m*counter1* sizeof(double));

    double * U   = (double *) malloc(n*counter*sizeof(double));
    //double * VT  = (double *) malloc(m*counter*sizeof(double));

    //transpose(V,VT, counter,m);
    //Transpose U
    transpose((double *)UT,U, n,counter);

#ifdef _MSC_VER
    _aligned_free(ATA);
    _aligned_free(AAT);
    _aligned_free(AT);
    _aligned_free(UT);
    //_aligned_free(V);
#else
    free(ATA);
    free(AAT);
    free(AT);
    free(UT);
    //free(S_2);
    //free(V);

#endif

    return new SVD{U,(double *)S,(double *)V,counter,n,m};

};

SVD * compute_svd(double * A,const  int m, const int n, double threshold, double accuracy){

#ifdef _MSC_VER
    const int N = 64;
    double * A_New = (double *) _aligned_malloc(n*m*sizeof(double));
#else
    const int N = 64;
    void * A_New;
    if(posix_memalign(&A_New , N, n*m      * sizeof(double))){
        cout<<"Error, couldn't allocate enough memory";
        exit(1);
    };
#endif

    memcpy(A_New,A,n*m* sizeof(double));
    SVD * svd =  power_methode_deflation((double *)A_New,m,n,accuracy, threshold);

#ifdef _MSC_VER
    _aligned_free(A_New);
#else
    free(A_New);
#endif

    return svd;
}

void free_svd_result(SVD* svd){
    delete svd;
}
