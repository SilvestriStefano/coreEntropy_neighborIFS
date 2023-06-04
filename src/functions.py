"""
Module that contains functions
"""
from sympy import Symbol, Function, evalf, Abs
from nested_lookup import nested_delete

import logging


# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s : %(message)s')
# add formatter to ch
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)


def check_vertices(
    edges:dict,
    vertices:dict,
    new_children:dict,
    current_h:str,
    current_h_val:complex,
    matching_h:list,
    nb_index:int,
    has_children:bool,
    label:str,
    R:complex,
    err:float,
    prec:int)->None:
    """
    It checks whether `current_h_val` is a new vertex and updates `edges` and `vertices` accordingly.
    It also updates `new_children` if necessary. 

    Parameters
    ----------
    edges: dict
        the dictionary of edges between the vertices
    vertices: dict
        the dictionary of vertices in the graph
    new_children: dict
        the dictionary of new child vertices in the graph
    current_h: str,
        the vertex being considered
    current_h_val: complex,
        the complex value of the vertex being considered
    matching_h: list,
        the list of the computed neighbors already present in the list of vertices
    nb_index: int,
        the current neighbor index
    has_children: bool,
        boolean value representing the existence of children of a vertex
    label: str,
        the edge label in the graph
    R: complex,
        the critical radius
    err: float,
        the allowed epsilon of error
    prec: int,
        the level of precision in evaluation

    Return
    ------
    None
    """
    
    if len(matching_h)==0: # current_h_val is POSSIBLY a new vertex
        if Abs(current_h_val).evalf(prec)<=R or Abs(Abs(current_h_val)-R).evalf(prec)<=err:
            has_children = False # current_h_val IS a child vertex
            nb_index+=1
            new_children.update({f"h{nb_index}": current_h_val})
            vertices.update({f"h{nb_index}": current_h_val})

            #if the current vertex has already some connections
            #update with a new one
            #otherwise create a new one
            if current_h in edges: 
                edges[current_h].update({f"h{nb_index}":{'label':label}})                        
            else:
                edges.update({current_h:{f"h{nb_index}":{'label':label}}})
        else: # current_h_val is NOT a VALID neighbor 
            has_children = True

    else: # current_h_val ALREADY EXISTS
        has_children = False 
        if current_h in edges:
            edges[current_h].update({matching_h[0]:{'label':label}})
        else:
            edges.update({current_h:{matching_h[0]:{'label':label}}})




def nbhG(param,maxDepth):
    """
    finds the edges in the neighbor graph for
    the parameter z.

    Parameters
    ----------
    param: complex number
        the parameter to check
    
    Attributes
    ----------
    maxDepth: int
        maximum depth
    
    Returns
    -------
    edges: dict
        the edges of the graph
    """
    z=Symbol('z')
    phiPM = Function('phiPM')(z)
    phiMP = Function('phiMP')(z)
    phiStar = Function('phiStar')(z)
    
    phiPM = (z-2)*param**(-1)# corresponds to fp^(-1) g fm
    phiMP = (z+2)*param**(-1)# corresponds to fm^(-1) g fp
    phiStar = z*param**(-1)# corresponds to fpm^(-1) g fpm
    
    # phiPM = (z-2*param)
    # phiMP = (z+2*param)
    # phiStar = (z)

    
    err = 1e-29
    prec = 30
    #initialize the dictionary of vertices in the graph
    vertices = {
        'id':0.,
        'h1':phiMP.evalf(prec,subs={z:0}) 
        # 'h1':phiPM.evalf(prec,subs={z:0}) #thanks to symmetry we can avoid this
    }
    #initialize the dictionary of new vertices at the current stage
    newVertices = {
        'h1':phiMP.evalf(prec,subs={z:0}) 
        # 'h1':phiPM.evalf(prec,subs={z:0}) #thanks to symmetry we can avoid this
    }
    #initialize the dictionary of edges between the vertices
    edges = {
        'id':{
              'h1':{'label':'- +','weight':0.25}
            # 'h1':{'label':'+ -','weight':0.75}#thanks to symmetry we can avoid this
            }
    }
    
    
    depth = 0
    neighborIndex = 1 #the label of the last vertex created
    
    criticalRad = (2*(1-Abs(param))**(-1)).evalf(prec) #the escape radius
    while len(newVertices) and depth<maxDepth:
        newChildren = {}
        verticesWithNoChild = {}
        
        #boolean values to check the existence of children of a vertex
        noChildPM = False
        noChildStar = False
        noChildMP = False
        
        logger.debug(f"{depth=} {newVertices=}")
        # ----------------------------- for each vertex in newVertices --------------------------------------------
        for keyNb,valNb in newVertices.items():
            logger.debug(f"{keyNb=} and {valNb=}")
            
            #compute the possible new neighbors
            hStar = phiStar.evalf(prec,subs={z:valNb})
            hPM = phiPM.evalf(prec,subs={z:valNb})
            hMP = phiMP.evalf(prec,subs={z:valNb})
            
            logger.debug(f"{hStar=}")
            logger.debug(f"{hPM=}")
            logger.debug(f"{hMP=}")

            #check if the computed neighbors exist already in the list of vertices
            matchStar = [key for key, value in vertices.items() if Abs(value-hStar).evalf(prec)<=err]
            matchPM = [key for key, value in vertices.items() if Abs(value-hPM).evalf(prec)<=err]
            matchMP = [key for key, value in vertices.items() if Abs(value-hMP).evalf(prec)<=err]

            # check_vertices(edges,vertices,newChildren,keyNb,hStar,matchStar,neighborIndex,noChildStar,' * ',criticalRad,err,prec)
            if len(matchStar)==0: # phiStar is POSSIBLY a new vertex
                logger.debug("phiStar is POSSIBLY a new vertex")
                if Abs(hStar).evalf(prec)<=criticalRad*(Abs(param)**(depth+2)) or Abs(Abs(hStar)-criticalRad*(Abs(param)**(depth+2))).evalf(prec)<=err:
                    logger.debug("phiStar IS a child vertex")
                    noChildStar = False # phiStar IS a child vertex
                    neighborIndex+=1
                    newChildren.update({f"h{neighborIndex}": hStar})
                    vertices.update({f"h{neighborIndex}": hStar})

                    #if the current vertex has already some connections
                    #update with a new one
                    #otherwise create a new one
                    if keyNb in edges: 
                        logger.debug(f"{keyNb=} is in edges, adding edge h{neighborIndex} with label *")
                        edges[keyNb].update({f"h{neighborIndex}":{'label':' * ', 'weight': 0.5}})                        
                    else:
                        logger.debug(f"{keyNb=} is NOT in edges, adding a NEW edge h{neighborIndex} with label *")
                        edges.update({keyNb:{f"h{neighborIndex}":{'label':' * ', 'weight': 0.5}}})
                else: # phiStar is NOT a VALID neighbor 
                    logger.debug("phiStar is NOT a new vertex")
                    noChildStar = True

            else: # phiStar ALREADY EXISTS
                logger.debug("phiStar ALREADY EXISTS")
                noChildStar = False 
                if keyNb in edges:
                    logger.debug(f"{keyNb=} is in edges, adding edge {matchStar[0]} with label *")
                    edges[keyNb].update({matchStar[0]:{'label':' * ', 'weight': 0.5}})
                else:
                    logger.debug(f"{keyNb=} is NOT in edges, adding a NEW edge {matchStar[0]} with label *")
                    edges.update({keyNb:{matchStar[0]:{'label':' * ', 'weight': 0.5}}})

            # check_vertices(edges,vertices,newChildren,keyNb,hPM,matchPM,neighborIndex,noChildPM,' + ',criticalRad,err,prec)
            if len(matchPM)==0: # phiPM is POSSIBLY a new vertex
                logger.debug("phiPM is POSSIBLY a new vertex")
                if Abs(hPM).evalf(prec)<=criticalRad*(Abs(param)**(depth+2)) or Abs(Abs(hPM)-criticalRad*(Abs(param)**(depth+2))).evalf(prec)<=err:
                    logger.debug("phiPM IS a child vertex")
                    noChildPM = False # phiPM IS a child vertex
                    neighborIndex += 1

                    newChildren.update({f"h{neighborIndex}": hPM})
                    vertices.update({f"h{neighborIndex}": hPM})

                    #if the current vertex has already some connections
                    #update with a new one
                    #otherwise create a new one
                    if keyNb in edges:
                        logger.debug(f"{keyNb=} is in edges, adding edge h{neighborIndex} with label + -")
                        edges[keyNb].update({f"h{neighborIndex}":{'label':'+ -', 'weight': 0.75}})
                    else:
                        logger.debug(f"{keyNb=} is NOT in edges, adding a NEW edge h{neighborIndex} with label + -")
                        edges.update({keyNb:{f"h{neighborIndex}":{'label':'+ -', 'weight': 0.75}}})
                else: # phiPM is NOT a VALID neighbor 
                    noChildPM = True

            else: # phiPM ALREADY EXISTS
                logger.debug("phiPM ALREADY EXISTS")
                noChildStar = False
                if keyNb in edges:
                    logger.debug(f"{keyNb=} is in edges, adding edge {matchPM[0]} with label + -")
                    edges[keyNb].update({matchPM[0]:{'label':'+ -', 'weight': 0.75}})
                else:
                    logger.debug(f"{keyNb=} is NOT in edges, adding a NEW edge {matchPM[0]} with label + -")
                    edges.update({keyNb:{matchPM[0]:{'label':'+ -', 'weight': 0.75}}})

            # check_vertices(edges,vertices,newChildren,keyNb,hMP,matchMP,neighborIndex,noChildMP,' - ',criticalRad,err,prec)
            if len(matchMP)==0: # phiMP is POSSIBLY a new vertex
                logger.debug("phiMP is POSSIBLY a new vertex")
                if Abs(hMP).evalf(prec)<=criticalRad*(Abs(param)**(depth+2)) or Abs(Abs(hMP)-criticalRad*(Abs(param)**(depth+2))).evalf(prec)<=err:
                    logger.debug("phiMP IS a child vertex")
                    noChildMP = False # phiMP IS a child vertex
                    neighborIndex += 1

                    newChildren.update({f"h{neighborIndex}": hMP})
                    vertices.update({f"h{neighborIndex}": hMP})

                    #if the current vertex has already some connections
                    #update with a new one
                    #otherwise create a new one
                    if keyNb in edges:
                        logger.debug(f"{keyNb=} is in edges, adding edge h{neighborIndex} with label - +")
                        edges[keyNb].update({f"h{neighborIndex}":{'label':'- +', 'weight': 0.25}})
                    else:
                        logger.debug(f"{keyNb=} is NOT in edges, adding a NEW edge h{neighborIndex} with label - +")
                        edges.update({keyNb:{f"h{neighborIndex}":{'label':'- +', 'weight': 0.25}}})
                else: # phiMP is NOT a VALID neighbor 
                    noChildMP = True

            else: # phiMP ALREADY EXISTS
                logger.debug("phiMP ALREADY EXISTS")
                noChildStar = False
                if keyNb in edges:
                    logger.debug(f"{keyNb=} is in edges, adding edge {matchMP[0]} with label - +")
                    edges[keyNb].update({matchMP[0]:{'label':'- +', 'weight': 0.25}})
                else:
                    logger.debug(f"{keyNb=} is NOT in edges, adding a NEW edge {matchMP[0]} with label - +")
                    edges.update({keyNb:{matchMP[0]:{'label':'- +', 'weight': 0.25}}})

            #in the case that all the computed neighbors are not valid
            #save the current neighbor
            if noChildStar and noChildPM and noChildMP:
                logger.debug(f"there are no new vertices. saving {keyNb=} in verticesWithNoChild")
                verticesWithNoChild.update({keyNb: valNb})
        #------------------------------------- end for loop ----------------------------------------------------
        
        #if there are neighbors without children
        #remove them from the list of vertices
        #and get rid of any edge connected to them 
        if len(verticesWithNoChild)!=0:
            vertices = {k:v for (k,v) in vertices.items() if k not in verticesWithNoChild }
            for key in verticesWithNoChild:
                edges = nested_delete(edges, key)
        
        #update the list of new vertices with the newly found vertices
        newVertices = newChildren
        depth += 1

    #remove the vertices that have no children
    #first find those with out-degree = 0 and remove them
    #i.e. those keys in `vertices`` that are not first-level key in `edges` 
    #then remove from `edges` the remaining first-level keys without properties
    #and all of its other instances
    nullOutDegre={k for k in vertices.keys() if k not in edges.keys()}
    for k in nullOutDegre:
        edges = nested_delete(edges,k)

    for k,v in edges.items():
        if len(v)==0:
            edges = nested_delete(edges,k)

    return edges


from src import angles as ang
from fractions import Fraction as Frac
from numpy import sort as npsort
from numpy import array as nparray
from numpy import ones as npones
from numpy import where as npwhere
from numpy import append as npappend

from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs

def core_entropy(*,num, den):
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
    logger.debug(csr_matrix((data, indices, indptr), shape=(len(tuples),len(tuples))).toarray())
    logger.debug(f"small evals using sparse {evals_small}")
    logger.debug(f"mid evals using sparse {evals_mid}")
    logger.debug(f"large evals using sparse {evals_large}")
    logger.debug(f"largest real eval = {max(set([x.real for x in [*evals_small,*evals_mid,*evals_large] if abs(x.imag)<.0001]))}")

    return max(set([x.real for x in [*evals_small,*evals_mid,*evals_large] if abs(x.imag)<.0001]))