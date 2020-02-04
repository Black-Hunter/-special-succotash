import os
# for easier development we can place the compilation command 
# in python itself (this compiles the module inplace)
os.system("cythonize -i partial_svd.pyx")

import numpy as np
from partial_svd import partial_svd

### testing implementation
# matrix for SVD
A = np.array([[1, 2, 3],
              [4, 5, 6],
              [7, 8, 9],
              [0, 11, 12]], dtype = np.double)
              

# now the partial svd
Ul, sl,VTl = partial_svd(A, 0.4)

print("\n")
print(Ul)
print('\n')
print(sl)
print('\n')
print(VTl)


# SVD from numpy
U, s, VT = np.linalg.svd(A, full_matrices=False)

#print(U,"\n", s,"\n", VT)
#print("\n")
#print(U)
#print('\n')
#print(s)
#print('\n')
#print(VT)
# test if correct
# atol -> absolute tolerance
#assert np.allclose(Ul @ np.diag(sl) @ VTl,
#                   U[:, 0:np.shape(Ul)[1]] @ np.diag(s[0:np.shape(sl)[0]]) @ VT[0:np.shape(VTl)[0], :], atol=1e-3)

