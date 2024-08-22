import numpy as np

def factor_mat_g(matG, nG, sB):
    """
    Reduces the positive definite matrix G to an upper triangular matrix.
    
    Parameters:
    matG: numpy.ndarray
        Semi-band of the matrix G, stored in a 2D numpy array.
    nG: int
        Order of the matrix G.
    sB: int
        Upper semi-band of the matrix G.
        
    Returns:
    numpy.ndarray
        The reduced matrix.
    """
    # Ensure matG is a numpy array
    matG = np.array(matG, dtype=np.complex128)
    
    # Iterate over columns of matG
    for j in range(nG):
        q = max(0, j - nG + sB)
        for i in range(sB - 1, q - 1, -1):
            k = j - i + sB
            if np.abs(matG[i, k]) != 0:
                c = matG[i, k] / matG[sB, j]
                m = sB
                n = k - 1
                for p in range(i, q - 1, -1):
                    m -= 1
                    n += 1
                    matG[m, n] -= c * matG[p, n]
                matG[i, k] = c
    
    return matG

# Example usage
# Initialize matrix and parameters
sB = 4
nG = 7
matG = np.zeros((sB, nG), dtype=np.complex128)

# Fill matG with example values
# (Replace this with actual values as needed)
matG[0, 2] = 13
matG[0, 3] = 24
matG[0, 4] = 35
matG[0, 5] = 46
matG[0, 6] = 57
matG[1, 1] = 12
matG[1, 2] = 23
matG[1, 3] = 34
matG[1, 4] = 45
matG[1, 5] = 56
matG[1, 6] = 67
matG[2, 0] = 11
matG[2, 1] = 22
matG[2, 2] = 33
matG[2, 3] = 44
matG[2, 4] = 55
matG[2, 5] = 66
matG[2, 6] = 77

# Perform factorization
reduced_matG = factor_mat_g(matG, nG, sB)
print(reduced_matG)

