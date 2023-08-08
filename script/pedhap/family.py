
import networkx as nx

class Family(nx.DiGraph):
    """ this is mostly a wrapper around networkx's DiGraph class
    
    Family members are stored as nodes, and children are linked to their parents
    by edges. More distant relatives must be identified by traversing the graph
    via child -> parent -> grandparent etc.
    """
    
    def __init__(self, id):
        self.id = id
        
        # properly initialise the networkx DiGraph class
        super(Family, self).__init__()
    
    def __repr__(self):
        return 'Family("{}")'.format(self.id)
    
    def __str__(self):
        return repr(self)
    
    def __gt__(self, other):
        return self.id > other.id
    
    def __iter__(self):
        # when we run though family members, only use members fully described
        # in the ped file (or inserted ourselves), otherwise we might get
        # individuals described only in child lines. This way, won't write ped
        # files with lines not present in the initial ped.
        return (x for x in iter(self._node) if not x._is_inferred() )
    
    def get_parents(self, person):
        return self.predecessors(person)
    
    def get_father(self, person):
        # return the male parent, or None if father not present
        for x in self.get_parents(person):
            if x.is_male():
                return x
    
    def get_mother(self, person):
        # return the female parent, or None if mother not present
        for x in self.get_parents(person):
            if not x.is_male() and not x.unknown_sex():
                return x
    
    def get_children(self, person):
        return self.successors(person)
    
    def add_person(self, person):
        if person is None:
            return
        if person.family != self.id:
            raise ValueError("{} didn't match family ID: {}".format(person.id,
                self.id))
        
        self.add_node(person)
    
    def set_mom(self, mom, child):
        if mom is None:
            return
        if mom not in self:
            raise ValueError("Can't set mother, not in family: {}".format(mom.id))
        
        if child not in self:
            raise ValueError("Can't set father, child not in family: {}".format(child.id))
        
        if self.get_mother(child) is not None and \
                self.get_mother(child) != mom:
            raise ValueError('adding a second mother to: {}'.format(child.id))
        
        if mom.id != child.mom:
            raise ValueError("mom ID not in child: {}".format(child.mom))
        
        # get the actual node for the mom, which contains the sex
        nodes = list(self.nodes)
        mom = nodes[nodes.index(mom)]
        if mom.is_male():
            raise ValueError("mom is not female: {}".format(mom.id))
        
        self.add_edge(mom, child)
    
    def set_dad(self, dad, child):
        if dad is None:
            return
        if dad not in self:
            raise ValueError("Can't set father, not in family: {}".format(dad.id))
            
        if child not in self:
            raise ValueError("Can't set father, child not in family: {}".format(child.id))
        
        if self.get_father(child) is not None and \
                self.get_father(child) != dad:
            raise ValueError('adding a second father to: {}'.format(child.id))
        
        if dad.id != child.dad:
            raise ValueError("dad ID not in child: {}".format(child.dad))
        
        # get the actual node for the dad, which contains the sex
        nodes = list(self.nodes)
        dad = nodes[nodes.index(dad)]
        if not (dad.is_male() or dad.unknown_sex()):
            raise ValueError("dad is not male: {}".format(dad.id))
        
        self.add_edge(dad, child)
    
    def get_proband(self):
        for x in self.nodes():
            if x.is_affected():
                return x
