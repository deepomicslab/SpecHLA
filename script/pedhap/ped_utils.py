
from pedigree import Trio
from person import Person
from typing import List
from family import Family

class Trio(object):
    def __init__(self, child, dad, mom) -> None:
        super().__init__()
        self.child = child
        self.dad = dad
        self.mom = mom
    def get_parent(self):
        return [self.dad, self.mom]
    def get_child(self):
        return self.child

    def __eq__(self, other):
        return (self.child.id == other.child.id and self.dad.id == other.dad.id and self.mom.id == other.mom.id)


def get_trios(family):
    """ get complete trios (child and parents) in a family
    
    Args:
        family: Family object (a graph with Persons as nodes)
    
    Returns:
        list of peds.Family objects, each for a unique trio
    """
    
    trios = []
    for x in family:
        mom = family.get_mother(x)
        dad = family.get_father(x)
        
        if mom is None and dad is None:
            continue
        
        # ignore people whose parents are based on their own ped line
        # if mom._is_inferred() or dad._is_inferred():
        #     continue
        
        trio = Trio(x, dad, mom)
        # trio.add_person(x)
        # trio.add_person(mom)
        # trio.add_person(dad)
        # trio.set_mom(mom, x)
        # trio.set_dad(dad, x)
        
        trios.append(trio)
    
    return trios

def get_top_level_trios(all_trios):
    r_trios = []
    for t in all_trios:
        flag = 1
        for i in t.get_parent():
            if i is not None and (i.mom != '0' or i.dad != '0'):
                flag = 0
                break
        if flag == 1:
            r_trios.append(t)
    return r_trios

def get_as_parent_trios(person: Person, all_trios: List[Trio]):
    r_trios = []
    for t in all_trios:
        for i in t.get_parent():
            if i.id == person.id:
                if t not in r_trios:
                    r_trios.append(t)
    return r_trios
def get_as_child_trio(person: Person, trios: List[Trio]):
    for t in trios:
        if t.get_child().id == person.id:
            return t

def get_next_level_trios(all_trios: List[Trio], current_level_trios: List[Trio]) -> List[Trio]:
    r_trios = []
    for t in current_level_trios:
        n_l_t = get_as_parent_trios(t.child, all_trios)
        for item in n_l_t:
            if item not in r_trios:
                r_trios.append(item)
    return r_trios

def get_bottom_level_trios(all_trios):
    r_trios = []
    for t in all_trios:
        if len(get_as_parent_trios(t.child, all_trios)) == 0:
            r_trios.append(t)
    return r_trios

def get_prev_level_trios(all_trios: List[Trio], current_level_trios: List[Trio]) -> List[Trio]:
    r_trios = []
    for t in current_level_trios:
        for i in t.get_parent():
            p_t = get_as_child_trio(i, all_trios)
            if p_t is not None:
                r_trios.append(p_t)
    return r_trios
def get_probands(family):
    """ find probands within a Family
    
    Returns:
        list of probands (as peds.Person objects)
    """
    probands = []
    for x in family:
        # arbitrarily define probands as individuals who are affected, and do
        # not have any children in the family. This could be a bit looser,
        # to cope with multigenerational familes, but will do for now.
        if x.is_affected():
            if len(list(family.get_children(x))) == 0:
                probands.append(x)
            else:
                mom = family.get_mom(x)
                dad = family.get_dad(x)
                
                if mom is None and dad is None:
                    continue
                
                # ignore people whose parents are based on their own ped line
                if mom._is_inferred() or dad._is_inferred():
                    continue
                
                probands.append(x)
    
    return probands
