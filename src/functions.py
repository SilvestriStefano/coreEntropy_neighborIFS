"""
Module that contains functions
"""
from sympy import Symbol, Function, evalf, Abs
from nested_lookup import nested_delete


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
    phiMP = (z+2)*param**(-1)# corresponds to fp^(-1) g fp
    phiStar = z*param**(-1)# corresponds to fpm^(-1) g fpm
    
    
    err = 1e-30
    prec = 30
    #initialize the dictionary of vertices in the graph
    vertices = {
        'id':0.,
        'h1':phiPM.evalf(prec,subs={z:0}),
        'h2':phiMP.evalf(prec,subs={z:0})
    }
    #initialize the dictionary of new vertices at the current stage
    newVertices = {
        'h1':phiPM.evalf(prec,subs={z:0}),
        'h2':phiMP.evalf(prec,subs={z:0})
    }
    #initialize the dictionary of edges between the vertices
    edges = {
        'id':{'h1':{'label':'+ -','weight':0.75},
              'h2':{'label':'- +','weight':0.25}}
    }
    
    
    depth = 0
    neighborIndex = 2; #the label of the last vertex created
    
    criticalRad = (2*(1-Abs(param))**(-1)).evalf(prec) #the escape radius
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
            hStar = phiStar.evalf(prec,subs={z:valNb})
            hPM = phiPM.evalf(prec,subs={z:valNb})
            hMP = phiMP.evalf(prec,subs={z:valNb})
            
            #check if the computed neighbors exist already in the list of vertices
            matchStar = [key for key, value in vertices.items() if Abs(value-hStar).evalf(prec)<=err]
            matchPM = [key for key, value in vertices.items() if Abs(value-hPM).evalf(prec)<=err]
            matchMP = [key for key, value in vertices.items() if Abs(value-hMP).evalf(prec)<=err]



            if len(matchStar)==0: # phiStar is POSSIBLY a new vertex
                if Abs(hStar).evalf(prec)<=criticalRad or Abs(Abs(hStar)-criticalRad).evalf(prec)<=err:
                    noChildStar = False # phiStar IS a child vertex
                    neighborIndex+=1
                    newChildren.update({f"h{neighborIndex}": hStar})
                    vertices.update({f"h{neighborIndex}": hStar})

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
                if Abs(hPM).evalf(prec)<=criticalRad or Abs(Abs(hPM)-criticalRad).evalf(prec)<=err:
                    noChildPM = False # phiPM IS a child vertex
                    neighborIndex += 1

                    newChildren.update({f"h{neighborIndex}": hPM})
                    vertices.update({f"h{neighborIndex}": hPM})

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
                if Abs(hMP).evalf(prec)<=criticalRad or Abs(Abs(hMP)-criticalRad).evalf(prec)<=err:
                    noChildMP = False # phiMP IS a child vertex
                    neighborIndex += 1

                    newChildren.update({f"h{neighborIndex}": hMP})
                    vertices.update({f"h{neighborIndex}": hMP})

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

    #remove the vertices that have no children
    for k,v in edges.items():
        if len(v)==0:
            edges = nested_delete(edges,k)

    return edges