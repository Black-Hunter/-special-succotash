#include "svd.h"
#include <iostream>
#include "omp.h"
#include "../headers/PerfEvent.hpp"


using namespace std;

int main(){
    int n = 3;
    int m = 3;

    double in [] = {12,22,33,44,55,66,77,88,99};
    SVD *svd;
    auto start = omp_get_wtime();

//    PerfEvent e;
//    {
//        e.startCounters();
        svd = compute_svd(in, n, m, 0.0, 0.00001);
//        e.stopCounters();
//        e.printReport(std::cout, n); // use n as scale factor
//    }

//    cout<<" Time Taken = "<< omp_get_wtime() - start<<"\n";


    int k = svd->k;
    m = svd->m;
    n = svd->n;
    const double * S = svd->s;
    const double *VT = svd->VT;
    const double * U = svd->U;

    cout<<" \n";

    for (int i = 0 ; i< k ; i++)
        cout<<S[i]<<" ";

    cout<<" \n";

    for (int i = 0 ; i< m*k ; i++)
        cout<<VT[i]<<" ";

    cout<<" \n";

    for (int i = 0 ; i< k*n ; i++)
        cout<<U[i]<<" ";


    free_svd_result(svd);

    return 0;

}
