from os import path
from typing import Union
from itertools import product
from src.neighbor import Neighbor
from sympy import Symbol, Function, Abs

import numpy as np
from numba import njit, prange
from numba.typed import List

import logging
import logging.config

src_dir, _ = path.split(path.abspath(__file__))
log_conf_path = path.join(path.dirname(src_dir),'log/logging.conf')
logging.config.fileConfig(log_conf_path)

# create logger
logger = logging.getLogger("default")

def is_child_neighbor(test_nbh:Neighbor, valid_nbhs:set, param:complex)->tuple:
    """
    Checks whether the Neighbor is a child Neighbor.
    
    Parameters
    ----------
    test_nbh: Neighbor
        the Neighbor to be tested
    valid_nbhs: set
        the set of valid Neighbors
    param: complex
        the parameter
    
    Return
    ------
    tuple of bool
    (is_new, is_child)
        (True, True) test_nbh is a new child Neighbor
        (True, False) test_nbh is not a Neighbor
        (False, True) test_nbh matches a Neighbor in the set
    """
    
    err = 1e-29
    prec = 30
    
    critical_rad = (2*(1-Abs(param))**(-1)).evalf(prec) #the escape radius
    is_new = test_nbh not in valid_nbhs
    
    if is_new: # test_nbh is POSSIBLY a new vertex
        logger.debug(f"{param:.5f};\t {test_nbh.word} is POSSIBLY a new neighbor")
        h_val = Abs(test_nbh.val)
        if h_val.evalf(prec)<=critical_rad or Abs(h_val-critical_rad).evalf(prec)<=err:
            logger.debug(f"{param:.5f};\t {test_nbh.word} IS a child vertex\n")
            is_child = True # test_nbh IS a child vertex
        else: # phi_Star is NOT a VALID neighbor 
            logger.debug(f"{param:.5f};\t {test_nbh.word} is NOT a new neighbor:\n\t\t {h_val=} {critical_rad=}\n")
            is_child = False
    else: # phi_Star ALREADY EXISTS
        logger.debug(f"{param:.5f}; \t {test_nbh.word} ALREADY EXISTS\n")
        is_child = True
    return (is_new,is_child)

def add_new_child(child_nbh:Neighbor, parent_nbh:Neighbor, edge:str, children:set, valid_nbhs:set, nbh_lookup:dict)->None:
    """Updates the set of child neighbors, the neighbor set and the dictionary of valid neighbors.
    
    Parameters
    ----------
    child_nbh: Neighbor
        the child Neighbor 
    parent_nbh: Neighbor
        the parent Neighbor 
    edge: str
        the function used to obtain child_nbh from parent_nbh
    valid_nbhs: set
        the set of valid Neighbors
    children: set
        The set of currently valid Neighbor children
    nbh_lookup: int
        The dictionary of the valid Neighbors. 
        It uses the Neighbor hash as the key and the Neighbor word as the value
    
    Return
    ------
    None
    """

    children.add(child_nbh)
    valid_nbhs.add(child_nbh)
    update_lookup(child_nbh,parent_nbh,edge,valid_nbhs,nbh_lookup,True)

def update_lookup(child_nbh:Neighbor, parent_nbh:Neighbor, edge:str, valid_nbhs:set, nbh_lookup:dict, is_new:bool)->None:
    """Updates the dictionary of valid neighbors and sets the relation between child and parent neighbor.
    
    Parameters
    ----------
    child_nbh: Neighbor
        the child Neighbor 
    parent_nbh: Neighbor
        the parent Neighbor 
    edge: str
        the function used to obtain child_nbh from parent_nbh
    valid_nbhs: set
        the set of valid neighbors
    nbh_lookup: int
        the dictionary of current neighbors
    is_new: bool
        whether child_nbh is a new valid neighbor
    
    Return
    ------
    None
    """ 
    if is_new:
        child_nbh_word = child_nbh.word
        nbh_lookup.update({child_nbh._hash:child_nbh_word})
    else:
        child_nbh_word = nbh_lookup[child_nbh._hash]
        
    valid_nbhs.remove(parent_nbh)
    parent_nbh.set_child(child_nbh_word,edge)
    valid_nbhs.add(parent_nbh)

def check_neighbor(test_nbh:Neighbor,
                   curr_nbh:Neighbor,
                   edge:str,
                   valid_nbhs:set,
                   children:set,
                   nbh_lookup:dict,
                   param:complex)->bool:
    """
    Paramteres
    ----------
    test_nbh: Neighbor
        The Neighbor to be checked
    curr_nbh: Neighbor
        The current Neighbor
    edge: str
        The string representation of the function that generates `test_nbh` from `curr_nbh`
    valid_nbh: set
        The set of currently valid Neighbors
    children: set
        The set of currently valid Neighbor children
    nbh_lookup: dict
        The dictionary of the valid Neighbors. 
        It uses the Neighbor hash as the key and the Neighbor word as the value
    param: complex
        The comple parameter that is being used.

    Return
    ------
    is_child : bool
    """
    is_new, is_child = is_child_neighbor(test_nbh,valid_nbhs,param)
    if is_new and is_child:
        add_new_child(test_nbh,curr_nbh,edge,children,valid_nbhs,nbh_lookup) 
    elif not is_new: 
        update_lookup(test_nbh,curr_nbh,edge,valid_nbhs,nbh_lookup,False)
    return is_child

def nbhG(param:complex, max_depth:int)->tuple:
    """
    Finds the edges in the Neighbor graph for
    the parameter z.

    Note that there might be Neighbors that are not valid.

    Parameters
    ----------
    param: complex
        the complex parameter to check
    maxDepth: int
        maximum depth
    
    Return
    ------
    valid_neighbors: set
        the set of valid Neighbors in the graph
    nbh_lookup: dict
        the dictionary of the valid Neighbors. 
        It uses the Neighbor hash as the key and the Neighbor word as the value
    """

    z = Symbol('z')
    phi_PM = Function('phiPM')(z)
    phi_MP = Function('phiMP')(z)
    phi_Star = Function('phiStar')(z)
    
    phi_PM = (z-2)*param**(-1)# corresponds to fp^(-1) g fm
    phi_MP = (z+2)*param**(-1)# corresponds to fm^(-1) g fp
    phi_Star = z*param**(-1)# corresponds to fpm^(-1) g fpm
        
    prec = 30
    
    #initialize the set of Neighbors in the graph
    valid_neighbors = set([
        Neighbor('.',0.+0.j,children=['+'],edges=['mp']),
        Neighbor('+',phi_MP.evalf(prec,subs={z:0}),parents=['.'])
    ])
    
    #initialize the dictionary of current Neighbors 
    nbh_lookup = {elem._hash:elem.word for elem in valid_neighbors}
    
    #initialize the set of new Neighbors at the current stage
    new_neighbors = set([Neighbor('+',phi_MP.evalf(prec,subs={z:0}),parents=['.'])])
    
    depth = 1
    
    while len(new_neighbors) and depth<max_depth:
        new_children = set()
        nbh_without_child = []
        
        #boolean values to check the existence of Neighbor's children
        is_child_Star = False
        is_child_PM = False
        is_child_MP = False
        
        for current_nbh in new_neighbors:
            current_word = current_nbh.word
            current_val = current_nbh.val
            
            #compute the possible new Neighbors
            h_Star = Neighbor(current_word+'0',phi_Star.evalf(prec,subs={z:current_val}),parents=[current_word])
            h_PM = Neighbor(current_word+'-',phi_PM.evalf(prec,subs={z:current_val}),parents=[current_word])
            h_MP = Neighbor(current_word+'+',phi_MP.evalf(prec,subs={z:current_val}),parents=[current_word])
            
            is_child_Star = check_neighbor(h_Star,current_nbh,'*',valid_neighbors,new_children,nbh_lookup,param)
            is_child_PM = check_neighbor(h_PM,current_nbh,'pm',valid_neighbors,new_children,nbh_lookup,param)
            is_child_MP = check_neighbor(h_MP,current_nbh,'mp',valid_neighbors,new_children,nbh_lookup,param)
                
            #in the case that all the computed neighbors are not valid
            #save the current Neighbor in a list 
            if not is_child_Star and not is_child_PM and not is_child_MP:
                logger.debug(f"{param:.5f}; {current_word} has no new child Neighbors")
                nbh_without_child.append(current_nbh)
            
        #if there are Neighbors without children
        #remove them from the set of valid Neighbors
        #and update the lookup dictionary
        if len(nbh_without_child)!=0:
            for elem in nbh_without_child:
                logger.debug(f"{param:.5f}; ...removing from valid Neighbors the ones with no children")
                valid_neighbors.remove(elem)
                valid_neighbors = {nbh.filter_children(elem.word) for nbh in valid_neighbors}
                del nbh_lookup[elem._hash]
        
                
        #update the list of new vertices with the newly found Neighbors
        new_neighbors.clear()
        new_neighbors.update(new_children)
        new_children.clear()
        depth += 1
        
    #clean up Neighbors with no children
    #NOTE: it might not find them all.
    logger.debug(f"{param:.5f}; ...another clean up of Neighbors without children")
    valid_neighbors = {nbh for nbh in valid_neighbors if len(nbh.children)>0}
    
    return valid_neighbors, nbh_lookup

def allsequences(n:int, terms:list=[1,-1] ,*,all:bool=False)->np.ndarray:
    """
    Generates from the elements in `terms` the list of all sequences 
    of length `n`.
    If `all=False` the sequences start with `terms[0]` and if 
    `terms=[1,0,-1]` it excludes the sequence `[1,0,0,0,..]`.

    Parameters
    ----------
    n: int
        Length of sequences
    terms: list
        Optional. Elements with which to construct the sequences.
        Default is `[1,-1]`.
    all: bool
        Decide wheter to include all possible sequences.
        Default is False. 

    Return
    ------
    numpy.array
        List of sequences.

    Example
    -------
    >>> seq_all_false = allsequences(3)
    >>> print(seq_all_f)
    array([[ 1, 1, 1],
           [ 1, 1,-1],
           [ 1,-1, 1],
           [ 1,-1,-1]])
    >>> seq_all_true = allsequences(3,all=True)
    >>> print(seq_all_true)
    array([[ 1,  1,  1],
           [ 1,  1, -1],
           [ 1, -1,  1],
           [ 1, -1, -1],
           [-1,  1,  1],
           [-1,  1, -1],
           [-1, -1,  1],
           [-1, -1, -1]])
    """
    
    if all:
        return np.array([seq for seq in product(terms, repeat=n)])
    
    lst = np.array([])
    for seq in product(terms, repeat=n):
        if np.all(np.equal(seq[1:],np.zeros(n-1))):
            continue
        if seq[0]==terms[0]: 
            lst = np.append(lst,seq)
        else:
            break
    return lst

@njit
def poly_eval(x:Union[int,float,complex],c:list)->Union[int,float,complex]:
    r"""Evaluate a polynomial at points x.
    If `c` is of length `n + 1`, this function returns the value
    .. math:: p(x) = c_0 + c_1 * x + ... + c_n * x^n
    
    Parameters
    ----------
    x: int, float, complex
        Value at which to evaluate the polynomial
    c: list
        List of coefficients of the polynomial
    
    Return
    ------
    int, float, complex
        Evaluation of the polynomial. The type depends on both `x` and `c`.
    """
    c0 = c[-1] + x*0
    for i in prange(2, len(c) + 1):
        c0 = c[-i] + c0*x
    return c0

@njit
def ps(pt:np.complex128, seqs:np.ndarray)->np.ndarray:
    """Evaluate at the given point the polynomial with coefficients 
    from the list of sequences.

    Parameters
    ----------
    pt: numpy.complex128
        Point at which to evaluate the power series
    seqs: numpy.ndarray
        Sequences of coefficients defining the polynomial
    
    Return
    ------
    numpy.ndarray
        Sequence of values
    """
    vals = []
    for i in prange(len(seqs)):
        vals.append(poly_eval(pt,seqs[i]))
    return np.asarray(vals)

@njit
def compute_green_MM0(pt_list:np.ndarray, level:int)->np.float64:
    r"""compute the minimum of the absolute vales from pt_list. 
    Take the log and normalize by the level.
    .. math:: \frac{1}{n} \log(\min (\left\vert \sum_{j=0}^{j=n-1}\epsilon_jx^j \right\vert))
    where 
    ..math:: \epsilon_j\in\{-1,1\}
    or 
    ..math:: \epsilon_j\in\{-1,0,1\}
    
    Parameters
    ----------
    pt_list: numpy.ndarray
        List of points
    level: int
        Level of approximation of the Green Function
    
    Return
    ------
    numpy.float64
        Value of the green function for the list of points.
    """
    vals = np.log(np.min(np.abs(pt_list)))/level
    return vals

@njit
def non_escaping_sequences(param:complex, sequences:list)->list:
    RAD = (1-np.abs(param))
    return [s for s in sequences if np.abs(1+param*poly_eval(param,s))*RAD < np.abs(param**(len(s)+1))]

def translate_str_to_int(i):
    CONDITIONS = [(lambda i: i=="+", 1), (lambda i: i=="-", -1), (lambda i: i=="0", 0)]
    for condition, replacement in CONDITIONS:
        if condition(i): return replacement
    return i

def get_coefficients(sequence:str)->list:
    return np.array(list(map(translate_str_to_int,sequence)))

@njit
def compare(sequence:list, s_list:list)->bool:
    n = len(sequence)
    for i,s in enumerate(s_list):
        s_trim = s[:n]
        if (s_trim==sequence).all():
            continue
        else:
            return False
    return True