
class Person(object):
    
    __slots__ = ('family', 'id', 'mom', 'dad', 'sex', 'phenotype', 'data', 'inferred')
    
    male_codes = set(['1', 'm', 'M', 'male'])
    female_codes = set(['2', 'f', 'F', 'female'])
    unknown_codes = set(['0', 'NA', 'unknown', '.', '-9'])
    
    def __init__(self, family, id, dad, mom, sex, phenotype, *args, **kwargs):
        
        self.family = family
        self.id = id
        self.mom = mom
        self.dad = dad
        self.sex = sex
        self.phenotype = phenotype
        self.data = args
        self.inferred = kwargs.get('inferred', False)
        # print(self.sex,"xxxxxx")
        
        if self.sex not in self.male_codes | self.female_codes | self.unknown_codes:
            raise ValueError('unknown sex code: {}'.format(self.sex))
        
        if self.phenotype not in set(['1', '2']) | self.unknown_codes:
            raise ValueError('unknown phenotype: {}'.format(self.phenotype))
    
    def __repr__(self):
        data = ''
        if len(self.data) > 0:
            temp = [ '"{}"'.format(x) for x in self.data ]
            data = ', {}'.format(", ".join(temp))
        
        return 'Person("{}", "{}", "{}", "{}", "{}", "{}"{})'.format(self.family,
            self.id, self.dad, self.mom, self.sex, self.phenotype, data)
    
    def __str__(self):
        """ convert the object back to a ped file line
        """
        
        data = ''
        if len(self.data) > 0:
            data = '\t' + '\t'.join(self.data)
        
        return '{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(self.family, self.id,
            self.dad, self.mom, self.sex, self.phenotype, data)
    
    def __hash__(self):
        return hash((self.family, self.id))
    
    def __eq__(self, other):
        return hash(self) == hash(other)
    
    def is_affected(self):
        return self.phenotype == "2"
    
    def is_male(self):
        return self.sex in self.male_codes
    
    def unknown_sex(self):
        return self.sex in self.unknown_codes

    def _is_inferred(self):
        return self.inferred
