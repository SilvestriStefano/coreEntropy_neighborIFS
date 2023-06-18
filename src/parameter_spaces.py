"""
Module that contains functions to generate the parameter spaces
Mandelbrot set
Thurston set
"""
import numpy as np
from numba import njit, prange
from src.utils import compute_green_MM0, ps, allsequences

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

def green_MM0(which:str, c:np.ndarray, level:int, x_dim:int = 500, y_dim:int = 500):
    """
    Computes the Green Function for the Barnsley or Thurston set according to
    the formula derived by Lindsey and Tiozzo in 2020.

    Parameters
    ----------
    which: str
        Choice of values 't' for Thurston set and 'b' for Barnsley set.
    c: numpy.ndarray
        array representing the complex plane
    level: int
        The level at which to approximate the Green's Function
    x_dim: int
        resolution for the real axis
    y_dim: int
        resolution for the imaginary axis

    Returns
    -------
    cntr: numpy.ndarray
        array of floats indicating the value of the Green Function at that particular parameter.
    """
    if which.lower() not in ['t','b']:
        raise ValueError("Only available options are `t` for Thurston set and `b` for Barnsley set")
    
    cntr = np.zeros_like(c, dtype=np.float64)
    seqs = allsequences(level,[1,0,-1]) if which.lower()=='b' else allsequences(level)
    for y in range(y_dim):
        for x in range(x_dim):
            if c[y][x]==0:
                cntr[y][x]=1e10
                continue
            if np.abs(c[y][x])>=2**(-0.25):
                cntr[y][x]=0
                continue
            pt = 1/c[y][x]
            cntr[y][x] = compute_green_MM0(ps(pt,seqs),level)
    return cntr

@np.vectorize
def green_MM0_contour(which:str, level:int, x:np.ndarray, y:np.ndarray)->np.ndarray:
    """
    Computes the Green Function for the Barnsley or Thurston set according to
    the formula derived by Lindsey and Tiozzo in 2020.

    It is meant to be used in combination with `matplotlib.pyplot.contour` and
    `matplotlib.pyplot.contourf` 

    Parameters
    ----------
    which: str
        Choice of values 't' for Thurston set and 'b' for Barnsley set.
    level: int
        The level at which to approximate the Green's Function
    x: numpy.ndarray
        The list of 'x' coordinate matrices
    y: numpy.ndarray
        The list of 'y' coordinate matrices

    Returns
    -------
    numpy.ndarray
        array of floats indicating the value of the Green Function at that particular parameter.
    
    Example
    -------
    >>> x_min, x_max,y_min, y_max = (0.001, 0.710,0.001, 0.710)
    >>> x_coords = np.linspace(x_min,x_max,x_dim)
    >>> y_coords = np.linspace(y_min,y_max,y_dim)
    >>> x, y = np.meshgrid(x_coords,y_coords)
    >>> z = green_MM0_contour(x,y)
    >>> plt.contourf(x, y, z)

    """
    if which.lower() not in ['t','b']:
        raise ValueError("Only available options are `t` for Thurston set and `b` for Barnsley set")
    seqs = allsequences(level,[1,0,-1]) if which.lower()=='b' else allsequences(level)
    return compute_green_MM0(ps(1/(x+y*1j),seqs),level)