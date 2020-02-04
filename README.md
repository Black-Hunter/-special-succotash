# SVD using the Power method
SVD Calculator using the Power method. Additionall for more time saving threshold  and accuracy could be added.

The Method can calculate not just all Egenvalues and Vectors but also all Eigenvalues which is bigger than the threshold. This means less calculation for just the needed Eigenvalues.  

Because the power Method iterate until it converges  to a specific Value, the ieteration can be also reduced. The number of iterations depends on the defference between the last calculated (inside the while loop --> see the method eigenvalue) Egenvalue and the current calculated. Usually it should be 0.0, but it also can be specifyed to make less iterations. This can be made by using the accuracy parameter.

## How to Use?
### Run The Programm:
To Use the SVD Method, you should at first include the header "svd.h" in your main. The Method to use in is _compute\_svd_ which delivers a SVD pointer and holds the two Matrixes V transposed and U, the vector S and the integers n,m and k.
The method  _compute\_svd_ takes five parameters, two of these parameters are optional namely the threshold and the accuracy. The following example demonstrate the Usage of The Method: 

```cpp
    
    // The input demensions
    int n = 3;
    int m = 3;
    double threshold = 0.0;
    double accuracy  = 0.00001;

    // The input Matrix
    double in [] = {
                    12,22,33,
                    44,55,66,
                    77,88,99
                   };
    
    // The results type 
    SVD *svd;

    // Computation
    svd = compute_svd(in, n, m, threshold, accuracy);

    // getting ther results
    int k = svd->k;
    m = svd->m;
    n = svd->n;
    const double * S = svd->s;
    const double *VT = svd->VT;
    const double * U = svd->U;
    
    // Memory cleaning
    free_svd_result(svd);

```

The perameter n from the code specify the rows and the parameter m defines the columens. The input Matrix can be from type one dimensional array or a pointer of type double. The method _free\_svd\_result_  is used to clean the results and free the Memory after getting the results.

### Test The Programm:    
For farthor implementation or modification the tests are written. They can be used to validate the new implementation. In the CMakeLists.txt the Targets are already defined.  

