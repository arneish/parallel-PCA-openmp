import numpy as NP
from scipy import linalg as LA
# A = [[108, 23, -33],[23, 121, -6],[-33, -6, 39]]
# A = NP.array(A, dtype=float)
# print(A)
# e_vals, e_v = LA.eig(A)
# print("e_vals:\n", e_vals)
# print("e_v:\n", e_v)

A = [[-1, 5, -9, 1],[9, 6, 0, 2],[2, -5, 1, 3]]
A = NP.array(A, dtype=float)
print(A)
u, s, vh = NP.linalg.svd(A, full_matrices=True)
print("u:\n", u)
print("sigma:\n", s)

temp = NP.matmul(u, s)
print(temp)
print("v_t:\n", vh)