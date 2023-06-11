from neighbor import Neighbor
from sympy import Symbol, Function, Abs

def is_child_neighbor(test_nbh:Neighbor, valid_nbhs:set, param:complex):
    """
    Checks whether the neighbor is a child neighbor.
    
    Parameters
    ----------
    test_nbh: Neighbor
        the Neighbor to be tested
    valid_nbhs: set
        the set of valid neighbors
    param: complex
        the parameter
    
    Returns
    -------
    tuple of bool
    (is_new, is_child)
        (True, True) test_nbh is a new child neighbor
        (True, False) test_nbh is not a neighbor
        (False, True) test_nbh matches a neighbor in the set
    """
    
    err = 1e-29
    prec = 30
    
    critical_rad = (2*(1-Abs(param))**(-1)).evalf(prec) #the escape radius
    is_new = test_nbh not in valid_nbhs
    for elem in valid_nbhs:
        diff = Abs(test_nbh.val-elem.val)
    
    if is_new: # test_nbh is POSSIBLY a new vertex
        # logger.debug(f"\t {test_nbh.word} is POSSIBLY a new neighbor")
        h_val = Abs(test_nbh.val)
        if h_val.evalf(prec)<=critical_rad or Abs(h_val-critical_rad).evalf(prec)<=err:
            # logger.debug(f"\t\t {test_nbh.word} IS a child vertex\n\n")
            is_child = True # test_nbh IS a child vertex
        else: # phi_Star is NOT a VALID neighbor 
            # logger.debug(f"\t\t {test_nbh.word} is NOT a new neighbor:\n\t\t {h_val=} {critical_rad=}\n\n")
            is_child = False
    else: # phi_Star ALREADY EXISTS
        # logger.debug(f"\t {test_nbh.word} ALREADY EXISTS\n\n")
        is_child = True
    return (is_new,is_child)

def add_new_child(child_nbh:Neighbor, parent_nbh:Neighbor, edge:str, children:set, valid_nbhs:set, nbh_lookup:dict):
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
        the set of valid neighbors
    children: list
        the set of valid neighbors
    nbh_lookup: int
        the dictionary of current neighbors
    """

    children.add(child_nbh)
    valid_nbhs.add(child_nbh)
    update_lookup(child_nbh,parent_nbh,edge,valid_nbhs,nbh_lookup,True)

def update_lookup(child_nbh:Neighbor, parent_nbh:Neighbor, edge:str, valid_nbhs:set, nbh_lookup:dict, is_new:bool):
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
    is_new, is_child = is_child_neighbor(test_nbh,valid_nbhs,param)
    if is_new and is_child:
        add_new_child(test_nbh,curr_nbh,edge,children,valid_nbhs,nbh_lookup) 
    elif not is_new: 
        update_lookup(test_nbh,curr_nbh,edge,valid_nbhs,nbh_lookup,False)
    return is_child

def nbhG(param:complex, max_depth:int)->dict:
    """
    Finds the edges in the neighbor graph for
    the parameter z.

    Note that there might be neighbors that are not valid.

    Parameters
    ----------
    param: complex
        the complex parameter to check
    maxDepth: int
        maximum depth
    
    Returns
    -------
    valid_neighbors: set
        the set of valid Neighbors in the graph
    nbh_lookup: dict
        the dictionary of the valid neighbors. 
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
    
    #initialize the set of neighbors in the graph
    valid_neighbors = set([
        Neighbor('.',0.+0.j,children=['+'],edges=['mp']),
        Neighbor('+',phi_MP.evalf(prec,subs={z:0}),parents=['.'])
    ])
    
    #initialize the dictionary of current neighbors 
    nbh_lookup = {elem._hash:elem.word for elem in valid_neighbors}
    
    #initialize the set of new neighbors at the current stage
    new_neighbors = set([Neighbor('+',phi_MP.evalf(prec,subs={z:0}),parents=['.'])])
    
    depth = 1
    
    while len(new_neighbors) and depth<max_depth:
        new_children = set()
        nbh_without_child = []
        
        #boolean values to check the existence of vertex's children
        is_child_Star = False
        is_child_PM = False
        is_child_MP = False
        
        for current_nbh in new_neighbors:
            current_word = current_nbh.word
            current_val = current_nbh.val
            
            #compute the possible new neighbors
            h_Star = Neighbor(current_word+'0',phi_Star.evalf(prec,subs={z:current_val}),parents=[current_word])
            h_PM = Neighbor(current_word+'-',phi_PM.evalf(prec,subs={z:current_val}),parents=[current_word])
            h_MP = Neighbor(current_word+'+',phi_MP.evalf(prec,subs={z:current_val}),parents=[current_word])
            
            is_child_Star = check_neighbor(h_Star,current_nbh,'*',valid_neighbors,new_children,nbh_lookup,param)
            is_child_PM = check_neighbor(h_PM,current_nbh,'pm',valid_neighbors,new_children,nbh_lookup,param)
            is_child_MP = check_neighbor(h_MP,current_nbh,'mp',valid_neighbors,new_children,nbh_lookup,param)
                
            #in the case that all the computed neighbors are not valid
            #save the current neighbor in a list 
            if not is_child_Star and not is_child_PM and not is_child_MP:
                nbh_without_child.append(current_nbh)
            
        #if there are neighbors without children
        #remove them from the set of valid neighbors
        #and update the lookup dictionary
        if len(nbh_without_child)!=0:
            for elem in nbh_without_child:
                valid_neighbors.remove(elem)
                valid_neighbors = {nbh.filter_children(elem.word) for nbh in valid_neighbors}
                del nbh_lookup[elem._hash]
        
                
        #update the list of new vertices with the newly found vertices
        new_neighbors.clear()
        new_neighbors.update(new_children)
        new_children.clear()
        depth += 1
        
    #clean up neighbors with no children
    #NOTE: it might not find them all.
    valid_neighbors = {nbh for nbh in valid_neighbors if len(nbh.children)>0}
    
    return valid_neighbors, nbh_lookup

# from nested_lookup import nested_delete
# def nbhG_old(param,maxDepth):
#     """
#     finds the edges in the neighbor graph for
#     the parameter z.

#     Parameters
#     ----------
#     param: complex number
#         the parameter to check
    
#     Attributes
#     ----------
#     maxDepth: int
#         maximum depth
    
#     Returns
#     -------
#     edges: dict
#         the edges of the graph
#     """
#     z=Symbol('z')
#     phiPM = Function('phiPM')(z)
#     phiMP = Function('phiMP')(z)
#     phiStar = Function('phiStar')(z)
    
#     phiPM = (z-2)*param**(-1)# corresponds to fp^(-1) g fm
#     phiMP = (z+2)*param**(-1)# corresponds to fm^(-1) g fp
#     phiStar = z*param**(-1)# corresponds to fpm^(-1) g fpm
    
#     # phiPM = (z-2*param)
#     # phiMP = (z+2*param)
#     # phiStar = (z)

    
#     err = 1e-29
#     prec = 30
#     #initialize the dictionary of vertices in the graph
#     vertices = {
#         'id':0.,
#         'h1':phiMP.evalf(prec,subs={z:0}) 
#         # 'h1':phiPM.evalf(prec,subs={z:0}) #thanks to symmetry we can avoid this
#     }
#     #initialize the dictionary of new vertices at the current stage
#     newVertices = {
#         'h1':phiMP.evalf(prec,subs={z:0}) 
#         # 'h1':phiPM.evalf(prec,subs={z:0}) #thanks to symmetry we can avoid this
#     }
#     #initialize the dictionary of edges between the vertices
#     edges = {
#         'id':{
#               'h1':{'label':'- +','weight':0.25}
#             # 'h1':{'label':'+ -','weight':0.75}#thanks to symmetry we can avoid this
#             }
#     }
    
    
#     depth = 0
#     neighborIndex = 1 #the label of the last vertex created
    
#     criticalRad = (2*(1-Abs(param))**(-1)).evalf(prec) #the escape radius
#     while len(newVertices) and depth<maxDepth:
#         newChildren = {}
#         verticesWithNoChild = {}
        
#         #boolean values to check the existence of children of a vertex
#         noChildPM = False
#         noChildStar = False
#         noChildMP = False
        
#         logger.debug(f"{depth=} {newVertices=}")
#         # ----------------------------- for each vertex in newVertices --------------------------------------------
#         for keyNb,valNb in newVertices.items():
#             logger.debug(f"{keyNb=} and {valNb=}")
            
#             #compute the possible new neighbors
#             hStar = phiStar.evalf(prec,subs={z:valNb})
#             hPM = phiPM.evalf(prec,subs={z:valNb})
#             hMP = phiMP.evalf(prec,subs={z:valNb})
            
#             logger.debug(f"{hStar=}")
#             logger.debug(f"{hPM=}")
#             logger.debug(f"{hMP=}")

#             #check if the computed neighbors exist already in the list of vertices
#             matchStar = [key for key, value in vertices.items() if Abs(value-hStar).evalf(prec)<=err]
#             matchPM = [key for key, value in vertices.items() if Abs(value-hPM).evalf(prec)<=err]
#             matchMP = [key for key, value in vertices.items() if Abs(value-hMP).evalf(prec)<=err]

#             # check_vertices(edges,vertices,newChildren,keyNb,hStar,matchStar,neighborIndex,noChildStar,' * ',criticalRad,err,prec)
#             if len(matchStar)==0: # phiStar is POSSIBLY a new vertex
#                 logger.debug("phiStar is POSSIBLY a new vertex")
#                 if Abs(hStar).evalf(prec)<=criticalRad*(Abs(param)**(depth+2)) or Abs(Abs(hStar)-criticalRad*(Abs(param)**(depth+2))).evalf(prec)<=err:
#                     logger.debug("phiStar IS a child vertex")
#                     noChildStar = False # phiStar IS a child vertex
#                     neighborIndex+=1
#                     newChildren.update({f"h{neighborIndex}": hStar})
#                     vertices.update({f"h{neighborIndex}": hStar})

#                     #if the current vertex has already some connections
#                     #update with a new one
#                     #otherwise create a new one
#                     if keyNb in edges: 
#                         logger.debug(f"{keyNb=} is in edges, adding edge h{neighborIndex} with label *")
#                         edges[keyNb].update({f"h{neighborIndex}":{'label':' * ', 'weight': 0.5}})                        
#                     else:
#                         logger.debug(f"{keyNb=} is NOT in edges, adding a NEW edge h{neighborIndex} with label *")
#                         edges.update({keyNb:{f"h{neighborIndex}":{'label':' * ', 'weight': 0.5}}})
#                 else: # phiStar is NOT a VALID neighbor 
#                     logger.debug("phiStar is NOT a new vertex")
#                     noChildStar = True

#             else: # phiStar ALREADY EXISTS
#                 logger.debug("phiStar ALREADY EXISTS")
#                 noChildStar = False 
#                 if keyNb in edges:
#                     logger.debug(f"{keyNb=} is in edges, adding edge {matchStar[0]} with label *")
#                     edges[keyNb].update({matchStar[0]:{'label':' * ', 'weight': 0.5}})
#                 else:
#                     logger.debug(f"{keyNb=} is NOT in edges, adding a NEW edge {matchStar[0]} with label *")
#                     edges.update({keyNb:{matchStar[0]:{'label':' * ', 'weight': 0.5}}})

#             # check_vertices(edges,vertices,newChildren,keyNb,hPM,matchPM,neighborIndex,noChildPM,' + ',criticalRad,err,prec)
#             if len(matchPM)==0: # phiPM is POSSIBLY a new vertex
#                 logger.debug("phiPM is POSSIBLY a new vertex")
#                 if Abs(hPM).evalf(prec)<=criticalRad*(Abs(param)**(depth+2)) or Abs(Abs(hPM)-criticalRad*(Abs(param)**(depth+2))).evalf(prec)<=err:
#                     logger.debug("phiPM IS a child vertex")
#                     noChildPM = False # phiPM IS a child vertex
#                     neighborIndex += 1

#                     newChildren.update({f"h{neighborIndex}": hPM})
#                     vertices.update({f"h{neighborIndex}": hPM})

#                     #if the current vertex has already some connections
#                     #update with a new one
#                     #otherwise create a new one
#                     if keyNb in edges:
#                         logger.debug(f"{keyNb=} is in edges, adding edge h{neighborIndex} with label + -")
#                         edges[keyNb].update({f"h{neighborIndex}":{'label':'+ -', 'weight': 0.75}})
#                     else:
#                         logger.debug(f"{keyNb=} is NOT in edges, adding a NEW edge h{neighborIndex} with label + -")
#                         edges.update({keyNb:{f"h{neighborIndex}":{'label':'+ -', 'weight': 0.75}}})
#                 else: # phiPM is NOT a VALID neighbor 
#                     noChildPM = True

#             else: # phiPM ALREADY EXISTS
#                 logger.debug("phiPM ALREADY EXISTS")
#                 noChildStar = False
#                 if keyNb in edges:
#                     logger.debug(f"{keyNb=} is in edges, adding edge {matchPM[0]} with label + -")
#                     edges[keyNb].update({matchPM[0]:{'label':'+ -', 'weight': 0.75}})
#                 else:
#                     logger.debug(f"{keyNb=} is NOT in edges, adding a NEW edge {matchPM[0]} with label + -")
#                     edges.update({keyNb:{matchPM[0]:{'label':'+ -', 'weight': 0.75}}})

#             # check_vertices(edges,vertices,newChildren,keyNb,hMP,matchMP,neighborIndex,noChildMP,' - ',criticalRad,err,prec)
#             if len(matchMP)==0: # phiMP is POSSIBLY a new vertex
#                 logger.debug("phiMP is POSSIBLY a new vertex")
#                 if Abs(hMP).evalf(prec)<=criticalRad*(Abs(param)**(depth+2)) or Abs(Abs(hMP)-criticalRad*(Abs(param)**(depth+2))).evalf(prec)<=err:
#                     logger.debug("phiMP IS a child vertex")
#                     noChildMP = False # phiMP IS a child vertex
#                     neighborIndex += 1

#                     newChildren.update({f"h{neighborIndex}": hMP})
#                     vertices.update({f"h{neighborIndex}": hMP})

#                     #if the current vertex has already some connections
#                     #update with a new one
#                     #otherwise create a new one
#                     if keyNb in edges:
#                         logger.debug(f"{keyNb=} is in edges, adding edge h{neighborIndex} with label - +")
#                         edges[keyNb].update({f"h{neighborIndex}":{'label':'- +', 'weight': 0.25}})
#                     else:
#                         logger.debug(f"{keyNb=} is NOT in edges, adding a NEW edge h{neighborIndex} with label - +")
#                         edges.update({keyNb:{f"h{neighborIndex}":{'label':'- +', 'weight': 0.25}}})
#                 else: # phiMP is NOT a VALID neighbor 
#                     noChildMP = True

#             else: # phiMP ALREADY EXISTS
#                 logger.debug("phiMP ALREADY EXISTS")
#                 noChildStar = False
#                 if keyNb in edges:
#                     logger.debug(f"{keyNb=} is in edges, adding edge {matchMP[0]} with label - +")
#                     edges[keyNb].update({matchMP[0]:{'label':'- +', 'weight': 0.25}})
#                 else:
#                     logger.debug(f"{keyNb=} is NOT in edges, adding a NEW edge {matchMP[0]} with label - +")
#                     edges.update({keyNb:{matchMP[0]:{'label':'- +', 'weight': 0.25}}})

#             #in the case that all the computed neighbors are not valid
#             #save the current neighbor
#             if noChildStar and noChildPM and noChildMP:
#                 logger.debug(f"there are no new vertices. saving {keyNb=} in verticesWithNoChild")
#                 verticesWithNoChild.update({keyNb: valNb})
#         #------------------------------------- end for loop ----------------------------------------------------
        
#         #if there are neighbors without children
#         #remove them from the list of vertices
#         #and get rid of any edge connected to them 
#         if len(verticesWithNoChild)!=0:
#             vertices = {k:v for (k,v) in vertices.items() if k not in verticesWithNoChild }
#             for key in verticesWithNoChild:
#                 edges = nested_delete(edges, key)
        
#         #update the list of new vertices with the newly found vertices
#         newVertices = newChildren
#         depth += 1

#     #remove the vertices that have no children
#     #first find those with out-degree = 0 and remove them
#     #i.e. those keys in `vertices`` that are not first-level key in `edges` 
#     #then remove from `edges` the remaining first-level keys without properties
#     #and all of its other instances
#     nullOutDegre={k for k in vertices.keys() if k not in edges.keys()}
#     for k in nullOutDegre:
#         edges = nested_delete(edges,k)

#     for k,v in edges.items():
#         if len(v)==0:
#             edges = nested_delete(edges,k)

#     return edges

