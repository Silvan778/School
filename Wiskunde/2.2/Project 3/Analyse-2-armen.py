import numpy as np
from scipy.linalg import eig
from numpy.linalg import matrix_rank

A = np.array([
 [ 0,  1, 0, 0, 0, 0],
 [0,  0, 0,  0, 0, 0],
 [0, 0, 0, 1, 0, 0],
 [0, 0, 21.42, 0, -21.34, 0],
 [ 0, 0, 0, 0, 0, 1],
 [ 0, 0, -28.58, 0, 77.62, 0]
])

B = np.array([
 [0],
 [ 0.09],
 [0],
 [-0.19],
 [0],
 [ 0.25]
])

def systeem_analyse(A, B):
    n = A.shape[0]
    
    print("Matrix A:")
    print(A)
    print("\nMatrix B:")
    print(B)
    
    # Stabiliteit
    eigenwaarden, eigenvectoren = eig(A)
    stabiel = np.all(np.real(eigenwaarden) < 0)
    print("\nEigenwaarden van A:", eigenwaarden)
    print("Stabiliteit:", stabiel)
    
    # Controleerbaarheid
    Ctrb = B
    for i in range(1, n):
        Ctrb = np.hstack((Ctrb, np.linalg.matrix_power(A, i) @ B))
    controleerbaar = matrix_rank(Ctrb) == n
    print("\nControleerbaarheidsmatrix:")
    print(Ctrb)
    print("Controleerbaarheid:", controleerbaar)
    
    # Stabiliseerbaarheid: focus op onstabiele modes
    onstabiel_indices = np.where(np.real(eigenwaarden) >= 0)[0]
    
    if len(onstabiel_indices) > 0:
        # Basis voor onstabiele modes
        onstabiele_eigenvectoren = eigenvectoren[:, onstabiel_indices]
        
        # Stabiliseerbaarheidsmatrix
        Stab_matrix = np.hstack([B, A @ B])
        for i in range(2, n):
            Stab_matrix = np.hstack([Stab_matrix, np.linalg.matrix_power(A, i) @ B])
        stabiliseerbaar = matrix_rank(Stab_matrix) == n
        print("\nStabiliseerbaarheidsmatrix:")
        print(Stab_matrix)
    else:
        Stab_matrix = np.array([])  # geen onstabiele modes
        stabiliseerbaar = True
        print("\nGeen onstabiele modes, stabiliseerbaarheidsmatrix leeg.")
    
    print("Stabiliseerbaarheid:", stabiliseerbaar)
    
    return {
        "Stabiliteit": stabiel,
        "Controleerbaarheid": controleerbaar,
        "Stabiliseerbaarheid": stabiliseerbaar
    }

resultaat = systeem_analyse(A, B)
print(resultaat)
