"""
Module that contains functions to generate the parameter spaces
Mandelbrot set
Thurston set
"""
import numpy as np
from numba import njit, prange

@njit
def mandelbrot(c, x_dim = 500, y_dim = 500, max_iter = 100):
    """
    It generates the Mandelbrot set.

    It is essentially the code from Murillo's group at MSU.

    Parameters
    ----------
    c: ndarray
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
    cntr = np.zeros_like(c, dtype=np.int64) #create a counter for each point in the plane
                                            #it will tell whether the parameter is in the mandelbrot set
    
    for y in prange(y_dim):
        for x in prange(x_dim):
            z = 0.0            #follow the orbit of the critical point
            it = 0             #count the iteration
            while (np.abs(z) <= max(2, np.abs(c[y][x])) ) and (it <= max_iter):
                z = z**2 + c[y][x]
                cntr[y][x] += 1
                it += 1
    return cntr

def checkPolyPZM(z, currVal, n, maxDeg):
    """
    Recursive function that checks whether the norm of any polynomial 
    with coefficients -1,0,+1 evaluated at z is larger than |z^n|/(1-|z|)
    
    Parameters
    ----------
    z: complex
        the input of the polynomial
    
    Arguments
    ---------
    currVal: int
        current value of the polynomial
    n: int
        current degree of the polynomial
    maxDeg: int
        maximum degree of the polynomial
    
    Returns
    -------
    result: int
        the degree of the polynomial
    """
    if n == maxDeg:
        return n
    
    #Check if this polynomial is too large
    if (np.abs(currVal) > (np.abs(z)**(n+1)/(1 - np.abs(z))) ):
        return n-1
    
    #Still good. Keep going: only keep the branch that doesn't die
    result = np.maximum(checkPolyPZM(z, currVal + z**(n+1), n+1, maxDeg), checkPolyPZM(z, currVal, n+1, maxDeg), checkPolyPZM(z, currVal - z**(n+1), n+1, maxDeg) )
    return result

def checkPolyPM(z, currVal, n, maxDeg):
    """
    Recursive function that checks whether the norm of any polynomial 
    with coefficients -1,0,+1 evaluated at z is larger than |z^n|/(1-|z|)
    
    Parameters
    ----------
    z: complex
        the input of the polynomial
    
    Arguments
    ---------
    currVal: int
        current value of the polynomial
    n: int
        current degree of the polynomial
    maxDeg: int
        maximum degree of the polynomial
    
    Returns
    -------
    result: int
        the degree of the polynomial
    """
    if n == maxDeg:
        return n
    
    #Check if this polynomial is too large
    if (np.abs(currVal) > (np.abs(z)**(n+1)/(1 - np.abs(z))) ):
        return n-1
    
    #Still good. Keep going: only keep the branch that doesn't die
    result = np.maximum(checkPolyPZM(z, currVal + z**(n+1), n+1, maxDeg), checkPolyPZM(z, currVal - z**(n+1), n+1, maxDeg) )
    return result