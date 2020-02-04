#include "catch.h"
#include "../Main/svd.h"
#include <iostream>
#include <memory>

using namespace std;
// Matrixes fpr Test
// not a Squared Matrix
const int n_1 = 4;
const int m_1 = 2;
double inM_1[] = {1,2,3,4,5,6,7,8};

// Squared Matrix
const int n_2 = 3;
const int m_2 = 3;
double   inM_2[] = {12,22,33,44,55,66,77,88,99};

// Validation Results
const int k_1 = 2;
const double  inM_1_T[] =  {1, 3, 5, 7, 2, 4, 6, 8};
const double  inM_AAT_1[] = {5, 11, 17, 23,11, 25, 39, 53,17, 39, 61, 83,23, 53, 83, 113};
const double   inM_1_S[] = {14.2274,1.257338};
const double   inM_1_V[] = {0.352062, 0.443626, 0.53519, 0.626754, 0.758981, 0.321242, 0.116498, 0.554238};
const double   inM_1_U[] = {0.376168, 0.926551 ,0.926551 ,0.376168};

// Squared Matrix
const int k_2 = 2;
const double   inM_2_T[] = {12, 44, 77, 22, 55, 88, 33, 66, 99};
const double   inM_AAT_2 []  = {1717, 3916, 6127, 3916, 9317, 14762, 6127, 14762, 23474};
const double   inM_2_S[] = {185.433, 11.0736, 0.176779};
const double   inM_2_V[] = {0.480462, 0.572035, 0.66478, 0.762758, 0.101553, 0.638661, 0.432847, 0.813918, 0.387532};
const double   inM_2_U[] = {0.217265, 0.874921, 0.432793, 0.520284, 0.271349, 0.809737, 0.825894, 0.401103, 0.396252};
const double   inM_2_U_sub[] = {0.217265, 0.825894, 0.271349 ,0.520284, 0.874921 ,0.401103};


bool isEqual(const double * m1, const double * m2, int size ){
    for(int i = 0; i < size; i++ )
    {
        if(abs(m1[i]) - abs(m2[i]) > 0.0001)
        {
            cout<<abs(m1[i]) << " m2 " << abs(m2[i])<<" m1 \n";
            return false;

        }
    }

    return true;

}


TEST_CASE("Matrix Operations", "[Matrix]"){

    double * T_1 = (double *) malloc(n_1*m_1*sizeof(double));
    double * T_2 = (double *) malloc(n_2*m_2*sizeof(double));

    double * AAT_1 = (double *) malloc(n_1*n_1*sizeof(double));
    double * AAT_2 = (double *) malloc(n_2*m_2*sizeof(double));

    transpose(inM_1, T_1, n_1,m_1 );
    transpose(inM_2, T_2, n_2,m_2 );

    matrix_matrixTranspose_Multi(inM_1, AAT_1, n_1, m_1, n_1);
    matrix_matrixTranspose_Multi(inM_2, AAT_2, n_2, m_2, n_2);

SECTION("Matrix Transpose")
{
REQUIRE(isEqual(inM_1_T, T_1,n_1*m_1) == 1);
    free(T_1);
REQUIRE(isEqual(inM_2_T, T_2,n_2*m_2) == 1);
    free(T_2);

}

SECTION("MATRIX MULTIPLICATION")
{
    REQUIRE(isEqual(AAT_1, inM_AAT_1,n_1*n_1) == 1);
    free(AAT_1);
    REQUIRE(isEqual(inM_AAT_2, AAT_2,n_2*n_2) == 1);
    free(AAT_2);

}



}


TEST_CASE("Results", "[SVD Result]"){
    // Quadratic Matrix
    SVD * result_1 = compute_svd(inM_1, 2, 4, -1);
    int k_1 = result_1->k;
    int m_1 = result_1->m;
    int n_1 = result_1->n;

    // None Quadratic Matrix
    SVD * result_2 = compute_svd(inM_2, 3,3, -1);
    int k_2 = result_2->k;
    int m_2 = result_2->m;
    int n_2 = result_2->n;

    // Threshhold = 10 --> k = 2
    SVD * result_3 = compute_svd(inM_2, 3,3, 10);
    int k_3 = result_3->k;
    int m_3 = result_3->m;
    int n_3 = result_3->n;

SECTION("VT Vector Result")
{
    REQUIRE(isEqual(result_1->VT,inM_1_V, k_1*n_1) == 1);
    REQUIRE(isEqual(result_2->VT,inM_2_V, k_2*n_2) == 1);
    REQUIRE(isEqual(result_3->VT,inM_2_V, k_3*n_3) == 1);

}

SECTION("U Matrix Result")
{
    REQUIRE(isEqual(result_1->U,inM_1_U, k_1*m_1) == 1);
    REQUIRE(isEqual(result_2->U,inM_2_U, k_2*m_2) == 1);
    REQUIRE(isEqual(result_3->U,inM_2_U_sub, k_3*m_3) == 1);

}

SECTION("S  Result")
{
    REQUIRE(isEqual(result_1->s,inM_1_S, k_1) == 1);
    REQUIRE(isEqual(result_2->s,inM_2_S, k_2) == 1);
    REQUIRE(isEqual(result_3->s,inM_2_S, k_3) == 1);

}

SECTION("Integers n,m,k")
{
    // N
    REQUIRE(result_1->n == 4);
    REQUIRE(result_2->n == 3);
    REQUIRE(result_3->n == 3);

    // M
    REQUIRE(result_1->m == 2);
    REQUIRE(result_2->m == 3);
    REQUIRE(result_3->m == 3);

    // N
    REQUIRE(result_1->k == 2);
    REQUIRE(result_2->k == 3);
    REQUIRE(result_3->k == 2);

}

}