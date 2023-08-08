
from person import Person
from family import Family

def open_ped(path):
    """ open a ped file and return a list of Family objects
    """
    
    # split file into lists of lines for each family, indexed by family ID
    fams = {}
    with open(path) as handle:
        sep = get_separator(handle)
        for line in handle:
            # ignore header and comment lines
            if line.startswith('family_id{}'.format(sep)) or line.startswith('#'):
                continue
            
            fam_id = line.split(sep, 1)[0]
            if fam_id not in fams:
                fams[fam_id] = Family(fam_id)
            
            person = Person(*line.strip().split(sep))
            if person in fams[fam_id]:
                raise ValueError('already family: {}'.format(person))
            
            fams[fam_id].add_person(person)
    
    return [ link_members(x) for x in list(fams.values()) ]

def get_separator(handle):
    """ get the column separator (assumes same on all lines)
    """
    current = handle.tell()
    line = handle.readline()
    handle.seek(current)
    
    if line.count('\t') > line.count(' '):
        return '\t'
    else:
        return ' '

def link_members(family):
    """ links family members, i.e. parents to children
    """
    
    # link parents to their children
    members = list(family)
    for person in members:
        # make placeholder parents, to match on family and individual IDs
        mom = Person(family.id, person.mom, 'NA', 'NA', 'female', 'NA', inferred=True)
        dad = Person(family.id, person.dad, 'NA', 'NA', 'male', 'NA', inferred=True)
        
        if mom.id != '0':
            if mom not in family:
                family.add_person(mom)
            family.set_mom(mom, person)
        
        if dad.id != '0':
            if dad not in family:
                family.add_person(dad)
            family.set_dad(dad, person)
    
    return family
