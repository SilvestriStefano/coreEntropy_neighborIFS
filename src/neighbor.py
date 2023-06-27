from sympy import Abs

class Neighbor:

    __slots__ = 'word','parents','children','val','edges','_hash'

    def __init__(self,
                 word:str,
                 val:complex,
                 *,
                 parents:list=None,
                 children:list=None,
                 edges:list=None):
        self.word = word
        self.val = complex(val)
        self.parents = parents if parents else []
        self.children = children if children else []
        self.edges = edges if edges else []
        self._hash = hash(complex(round(self.val.real,13),round(self.val.imag,13)))

    def __repr__(self)->str:
        return f"Neighbor({self.word})"
    
    def __hash__(self)->int:
        return self._hash
    
    def __eq__(self, other)->bool:
        if not isinstance(other, type(self)): return NotImplemented
        return Abs(self.val-other.val)<=1e-13
    
    def set_parent(self, elem:str)->None:
        """Set the Neighbor's parent
        
        Parameters
        ----------
        elem: str
          The itinerary of the parent Neighbor
        
        Returns
        -------
        None
        """
        self.parents.append(elem)
    
    def set_child(self, elem:str, edge:str)->None:
        """Set the Neighbor's child
        
        Parameters
        ----------
        elem: str
          The itinerary of the child Neighbor
        edge: str
          The string representing how to obtain the child Neighbor
        
        Returns
        -------
        None
        """
        self.children.append(elem)
        self.edges.append(edge)
    
    def filter_children(self, child_to_filter:str):
        """Remove a specified child from the Neighbor's children
        
        Parameters
        ----------
        child_to_filter: str
          The itinerary of the child Neighbor to be removed
        
        Returns
        -------
        The Neighbor's instance
        """
        self.children=[child for child in self.children if child != child_to_filter]
        return self