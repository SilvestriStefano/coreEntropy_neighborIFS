"""
Module that contains functions to generate the dynamical spaces
Julia set
Attractor of an IFS
"""
import numpy as np
from numba import njit, prange

@njit 
def julia(c, z, x_dim = 500, y_dim = 500, max_iter = 100):
    """
    It generates the Julia set.

    It is essentially the code from Murillo's group at MSU.

    Parameters
    ----------
    c: complex
        parameter defining the quadratic map
    z: ndarray
        array representing the complex plane

    Arguments
    ---------
    x_dim: int
        resolution for the real axis
    y_dim: int
        resolution for the imaginary axis
    max_iter: int
        maximum number of iteration

    Returns
    -------
    cntr: ndarray
        array of integers indicating how long it took for 0 to escape for that particular parameter
    """
    cntr = np.zeros_like(z, dtype=np.int64) #create a counter for each point in the plane
                                            #it will tell whether the parameter is in the filled julia set 
    
    for y in prange(y_dim):
        for x in prange(x_dim):
            it = 0
            while (np.abs(z[y][x]) <= 2) and (it <= max_iter):
                z[y][x] = z[y][x]**2 + c #iterate 
                cntr[y][x] += 1                      
                it += 1                          
    return cntr