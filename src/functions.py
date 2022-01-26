"""
Module that contains functions
"""
import numpy as np
from mpmath import mp
from nested_lookup import nested_delete


def neighGraphAlg(z,maxDepth):
    """
    finds the edges in the neighbor graph for
    the parameter z.

    Parameters
    ----------
    z: complex number
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
    prec = 1e-14
    #initialize the dictionary of vertices in the graph
    vertices = {
        'id':0.,
        'h1':mp.mpc(-2./param),
        'h2':mp.mpc(2./param)
    }
    #initialize the dictionary of new vertices at the current stage
    newVertices = {
        'h1':mp.mpc(-2./param),
        'h2':mp.mpc(2./param)
    }
    #initialize the dictionary of edges between the vertices
    edges = {
        'id':{'h1':{'label':'+ -','weight':0.75},
              'h2':{'label':'- +','weight':0.25}}
    }
    
    
    depth = 0
    neighborIndex = 2; #the label of the last vertex created
    
    criticalRad = mp.mpf(2./(1-np.abs(param))) #the escape radius
    
    while len(newVertices) and depth<maxDepth:
        newChildren = {}
        verticesWithNoChild = {}
        
        #boolean values to check the existence of children of a vertex
        noChildPM = False
        noChildStar = False
        noChildMP = False
        
        # ----------------------------- for each vertex in newVertices --------------------------------------------
        for keyNb,valNb in newVertices.items():
            #compute the possible new neighbors
            phiStar = mp.mpc(1/param*(valNb)) # corresponds to fpm^(-1) g fpm
            phiPM = mp.mpc(1/param*(valNb-2)) # corresponds to fp^(-1) g fm
            phiMP = mp.mpc(1/param*(valNb+2)) # corresponds to fm^(-1) g fp
            
            #check if the computed neighbors exist already in the list of vertices
            matchStar = [key for key, value in vertices.items() if np.abs(value-phiStar) <= prec]#if value==phiStar] #
            matchPM = [key for key, value in vertices.items() if np.abs(value-phiPM) <= prec]#if value==phiPM] #
            matchMP = [key for key, value in vertices.items() if np.abs(value-phiMP) <= prec]#if value==phiMP] #
            
            
            if len(matchStar)==0: # phiStar is POSSIBLY a new vertex
                if mp.mpf(np.abs(phiStar))<=criticalRad or mp.mpf(np.abs(np.abs(phiStar)-criticalRad))<=prec:
                    noChildStar = False # phiStar IS a child vertex
                    neighborIndex+=1
                    
                    newChildren.update({f"h{neighborIndex}": phiStar})
                    vertices.update({f"h{neighborIndex}": phiStar})
                    
                    #if the current vertex has already some connections
                    #update with a new one
                    #otherwise create a new one
                    if keyNb in edges: 
                        edges[keyNb].update({f"h{neighborIndex}":{'label':' * ', 'weight': 0.5}})                        
                    else:
                        edges.update({keyNb:{f"h{neighborIndex}":{'label':' * ', 'weight': 0.5}}})
                else: # phiStar is NOT a VALID neighbor 
                    noChildStar = True
                    
            else: # phiStar ALREADY EXISTS
                noChildStar = False 
                if keyNb in edges:
                    edges[keyNb].update({matchStar[0]:{'label':' * ', 'weight': 0.5}})
                else:
                    edges.update({keyNb:{matchStar[0]:{'label':' * ', 'weight': 0.5}}})
                
            if len(matchPM)==0: # phiPM is POSSIBLY a new vertex
                if mp.mpf(np.abs(phiPM))<=criticalRad or mp.mpf(np.abs(np.abs(phiPM)-criticalRad))<=prec:
                    noChildPM = False # phiPM IS a child vertex
                    neighborIndex += 1
                                        
                    newChildren.update({f"h{neighborIndex}": phiPM})
                    vertices.update({f"h{neighborIndex}": phiPM})
                    
                    #if the current vertex has already some connections
                    #update with a new one
                    #otherwise create a new one
                    if keyNb in edges:
                        edges[keyNb].update({f"h{neighborIndex}":{'label':'+ -', 'weight': 0.75}})
                    else:
                        edges.update({keyNb:{f"h{neighborIndex}":{'label':'+ -', 'weight': 0.75}}})
                else: # phiPM is NOT a VALID neighbor 
                    noChildPM = True
            
            else: # phiPM ALREADY EXISTS
                noChildStar = False
                if keyNb in edges:
                    edges[keyNb].update({matchPM[0]:{'label':'+ -', 'weight': 0.75}})
                else:
                    edges.update({keyNb:{matchPM[0]:{'label':'+ -', 'weight': 0.75}}})
            
            if len(matchMP)==0: # phiMP is POSSIBLY a new vertex
                if mp.mpf(np.abs(phiMP))<=criticalRad or mp.mpf(np.abs(np.abs(phiMP)-criticalRad))<=prec:
                    noChildMP = False # phiMP IS a child vertex
                    neighborIndex += 1
                    
                    newChildren.update({f"h{neighborIndex}": phiMP})
                    vertices.update({f"h{neighborIndex}": phiMP})
                    
                    #if the current vertex has already some connections
                    #update with a new one
                    #otherwise create a new one
                    if keyNb in edges:
                        edges[keyNb].update({f"h{neighborIndex}":{'label':'- +', 'weight': 0.25}})
                    else:
                        edges.update({keyNb:{f"h{neighborIndex}":{'label':'- +', 'weight': 0.25}}})
                else: # phiMP is NOT a VALID neighbor 
                    noChildMP = True
            
            else: # phiMP ALREADY EXISTS
                noChildStar = False
                if keyNb in edges:
                    edges[keyNb].update({matchMP[0]:{'label':'- +', 'weight': 0.25}})
                else:
                    edges.update({keyNb:{matchMP[0]:{'label':'- +', 'weight': 0.25}}})
                    
            #in the case that all the computed neighbors are  not valid
            #save the current neighbor
            if noChildStar and noChildPM and noChildMP:
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
        
    
    return edges