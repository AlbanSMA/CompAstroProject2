import numpy as np

#rate of the different transitions
r_1 = 1e-43
r_2 = 1e-26
r_3 = 1e-34

#densities
n_p = 1e26
n_2He = 1e26
n_3He = 1e26
n_4He = 0.01e26

#values within the matrix
J_11 = 2*r_1*n_p + r_2 * n_2He
J_12 = r_2*n_p
J_13 = -2*r_3*n_3He
J_14 = 0

J_21 = -2*r_1*n_p + r_2*n_2He
J_22 = r_2*n_p
J_23 = 0
J_24 = 0

J_31 = -r_2*n_2He
J_32 = -r_2*n_p
J_33 = 2*r_3*n_3He
J_34 = 0

J_41 = 0
J_42 = 0
J_43 = -2*r_3*n_3He
J_44 = 0

J = np.array([J_11, J_12, J_13, J_14, J_21, J_22, J_23, J_24,
              J_31, J_32, J_33, J_34, J_41, J_42, J_43, J_44]).reshape(4,4)

eigenvalues, eigenvectors = np.linalg.eig(J)
print(eigenvalues)

stiffness = max(eigenvalues)/min(eigenvalues)
print(stiffness)

