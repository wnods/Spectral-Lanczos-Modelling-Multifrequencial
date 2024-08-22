import numpy as np

def solve_ma_gx(matG, vecG, nG, sB):
    """
    Reduces the source vector vecG and calculates the solution of the system (back substitution).
    
    Parameters:
    matG: numpy.ndarray
        Semi-band of the global matrix G, stored in a 2D numpy array.
    vecG: numpy.ndarray
        Source vector to be solved for.
    nG: int
        Order of the matrix G.
    sB: int
        Upper semi-band of the matrix G.
        
    Returns:
    numpy.ndarray
        The solution vector. Note that the source vector is replaced by the solution.
    """
    # Ensure matG and vecG are numpy arrays of complex numbers
    matG = np.array(matG, dtype=np.complex128)
    vecG = np.array(vecG, dtype=np.complex128)
    
    # Forward substitution
    for j in range(nG):
        k = max(0, j - nG + sB)
        for m in range(sB - 1, k - 1, -1):
            n = j - m + sB
            vecG[n] -= matG[m, n] * vecG[j]
        vecG[j] /= matG[sB - 1, j]  # Adjust index for 0-based Python
    
    # Back substitution
    for j in range(1, nG):
        i = nG - 1 - j
        k = max(0, i - nG + sB)
        for m in range(sB - 1, k - 1, -1):
            n = i - m + sB
            vecG[i] -= matG[m, n] * vecG[n]
    
    return vecG

# Example usage
sB = 4
nG = 7
matG = np.zeros((sB, nG), dtype=np.complex128)
vecG = np.zeros(nG, dtype=np.complex128)

# Fill matG and vecG with example values
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

# Example source vector (replace with actual values)
vecG = np.array([1+0j, 2+0j, 3+0j, 4+0j, 5+0j, 6+0j, 7+0j], dtype=np.complex128)

# Solve the system
solution_vecG = solve_ma_gx(matG, vecG, nG, sB)
print(solution_vecG)

