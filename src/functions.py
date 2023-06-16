"""
Module that contains functions
"""
from typing import Union
from pathlib import Path
from src.utils import nbhG

from src import angles as ang
from fractions import Fraction as Frac
from numpy import array as nparray
from numpy import ones as npones
from numpy import where as npwhere
from numpy import append as npappend
from numpy import float64

from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs

import logging
import logging.config

log_conf_path = Path('.') / 'logging.conf'
logging.config.fileConfig(log_conf_path)

# create logger
logger = logging.getLogger("default")


def create_nbh_graph(param:Union[int,float,complex], max_depth:int)->dict:
    """
    Creates the Neighbor Graph following NetworkX's graph data structure:
    a 'dictionary of dictionaries of dictionaries'. 
    
    For param values with absolute value less than 0.5, it returns an empty dictionary.
    For param values with absolute value bigger than 0.5 but whose associated
    attractor is a Cantor set, it might return a non-empty dictionary. However there should
    not be any loop.
    For param values whose associate attractor is connected, there should be at least one loop.

    Note
    ----
    It might contain non valid Neighbors, i.e. vertices with no children.

    Parameters
    ----------
    param: int, float, complex
      A number with absolute value less than 1.
    max_depth: int
      How long should the algorithm continue before exiting.
    
    Returns
    -------
    dict
      The Neighbor Graph

    Example
    -------
    >>> G = create_nbh_graph(0.5+0.*1j,6)
    >>> print(G)
    {'id': {'h+': 'mp'}, 'h+':{'h+': 'pm'}}
    """
    try:
        param = complex(param)
    except ValueError:
        raise ValueError("The parameter should be a complex number")
    
    if type(max_depth) is not int:
        raise ValueError("The maximum depth should be an integer value")
    
    if abs(param)<0.5 or abs(abs(param)-1.0)<1e-13 or abs(param)>1:
        return {}

    valid_nbh, nbh_lookup = nbhG(param,max_depth)
    nbh_graph = {}
    for nbh in valid_nbh:
        if nbh.word==".": 
            vertex_label="id"
        else:
            vertex_label = f"h{nbh_lookup[nbh._hash]}"
        connected_to = {f"h{child}":edge for child,edge in zip(nbh.children,nbh.edges)}
        nbh_graph.update({vertex_label:connected_to})
    return nbh_graph


def core_entropy(*,num:int, den:int)->float64:
    """
    Calculates the core entropy for a given rational angle. 
    
    Parameters
    ----------
    num: int
      The numerator of the rational angle.
    den: int
      The denominator of the rational angle.
    
    Returns
    -------
    numpy.float64
      The core entropy
    """
    
    if type(num) is not int or type(den) is not int:
        raise ValueError("Arguments should be integers or strings of integer")
    if type(int(num)) is not int or type(int(den)) is not int:
        raise ValueError("Arguments should be integers or strings of integer")

    #define a new rational angle
    theta = ang.Angle(num,den)
    
    thetaFr = theta.frac #represent it as a fraction

    if thetaFr==Frac(1,2):
        return 2.0
    if thetaFr==Frac(0,1) or thetaFr==Frac(1,1) :
        return 1.0
   
    #compute the orbit and find the period
    period_length, period_start = theta.period()
    orb = theta.orbit_list
    preperiod_length = len(orb)-1-period_length
    
    #partition of the circle
    intOne = (thetaFr*Frac(1,2),(thetaFr+Frac(1,1))*Frac(1,2))
    intTwo = (intOne[1],intOne[0])
    
    logger.debug(f"The angle is {thetaFr} \n which has a preperiod of {preperiod_length} and a period of {period_length} \n and the orbit is {orb}\n")
    logger.debug(f"partition the circle in two intervals: {intOne} and {intTwo}\n") #{(intOne[1],intOne[0])}\n") #
    
    #create the vertex set of the wedge
    tuples = nparray(["%d-%d"%(i,j) for i in range(1,len(orb)) for j in range(1,len(orb)) if i<j])
    logger.debug(f"tuples is {tuples}\n")
    
    #define the separated and non-separated vertices
    dicTuples = {}
    entries = 0
    indices = nparray([],dtype=int)
    indptr = nparray([0],dtype=int)
    
    for tup in tuples:
        i,j=tup.split('-')
        i=int(i)
        j=int(j)
        max_ind = len(orb)-1
        
        logger.debug(f"the tuple {tup} corresponding to the angles {orb[i-1]} and {orb[j-1]}")

        if (
            (intOne[0]<=orb[i-1]<intOne[1] and intOne[0]<=orb[j-1]<intOne[1]) or 
            ((intTwo[0]<=orb[i-1] or orb[i-1]<intTwo[1]) and (intTwo[0]<=orb[j-1] or orb[j-1]<intTwo[1]))):
            logger.debug(f"the tuple {tup} is not separated.")
            
            if (j+1)<=max_ind:
                target = str(i+1)+'-'+str(j+1)
            else:
                new_j = period_start+1 
                new_i = i+1
                target = str(new_i)+'-'+str(new_j) if new_i<new_j else str(new_j)+'-'+str(new_i)
            logger.debug(f"target = {target}")
            try:
                index = npwhere(tuples==target)[0][0]
            except IndexError:
                logger.error(f"there is no tuple {target}")
            else:
                logger.debug(f"the target is {target} which has index {index}\n")
                indices = npappend(indices,index)

            dicTuples.update({tup:{'sep':False,'mapsTo':[target]}})
            entries += 1
        else:
            logger.debug(f"the tuple {tup} is separated")
            target_one = str(1)+'-'+str(i+1) if (i+1)<=max_ind else str(1)+'-'+str(period_start+1 if (period_start+1)!=1 else 0 )
            logger.debug(f"target one = {target_one}")

            target_two = str(1)+'-'+str(j+1) if (j+1)<=max_ind else str(1)+'-'+str(period_start+1 if (period_start+1)!=1 else 0 )
            logger.debug(f"target two = {target_two}")
            
            try: 
                index_one = npwhere(tuples==target_one)[0][0]
            except IndexError:
                logger.error(f"there is no tuple {target_one}")
            else:
                logger.debug(f"the target_one is {target_one} which has index {index_one}")
                indices = npappend(indices,index_one)
                entries+=1
            try:
                index_two = npwhere(tuples==target_two)[0][0]
            except IndexError:
                logger.error(f"there is no tuple {target_two}")
            else:
                logger.debug(f"the target_two is {target_two} which has index {index_two}\n")
                indices = npappend(indices,index_two)
                entries+=1
            
            
            dicTuples.update({tup:{'sep':True,'mapsTo':[target_one,target_two]}})
            # entries += 2
        indptr = npappend(indptr,entries)
    
    data = npones(entries,dtype='float64')
    
    logger.debug(dicTuples)
    logger.debug(f"{data=}")
    logger.debug(f"{indices=}")
    logger.debug(f"{indptr=}")
    
    kE = min(2,len(tuples)-2)
    
    adj_matrix = csr_matrix((data, indices, indptr), shape=(len(tuples),len(tuples)))
    try:
        # evals_small = eigs(adj_matrix,k=kE, sigma=0.00000001, which='LM',return_eigenvectors=False)
        # evals_mid = eigs(adj_matrix,k=kE, sigma=0.8999999, which='LM',return_eigenvectors=False)
        evals_large = eigs(adj_matrix,k=kE, sigma=1.7999999, which='LM',return_eigenvectors=False)
    except:
        return 1.0
    else:
        return max(set([x.real for x in evals_large if abs(x.imag)<.0001]))
    # logger.debug(csr_matrix((data, indices, indptr), shape=(len(tuples),len(tuples))).toarray())
    # logger.debug(f"small evals using sparse {evals_small}")
    # logger.debug(f"mid evals using sparse {evals_mid}")
    # logger.debug(f"large evals using sparse {evals_large}")
    # logger.debug(f"largest real eval = {max(set([x.real for x in [*evals_small,*evals_mid,*evals_large] if abs(x.imag)<.0001]))}")

    # return max(set([x.real for x in [*evals_small,*evals_mid,*evals_large] if abs(x.imag)<.0001]))