'''
Created on Nov 19, 2016

@author: Lucy
'''

from gurobipy.gurobipy import tupledict

class IncidenceMatrix(object):
    '''
    classdocs
    '''

    def __init__(self, reactions):
        '''
        reactions: A list of cobra reactions
        metabolites: A list of the cobra metabolites
        '''
        self._reactions = reactions
        self._metabolites = []
        
        self._mat = tupledict()
        for r in reactions:
            for m in r.metabolites:
                self._mat[r.id, m.id] = r.get_coefficient(m)
                if m not in self._metabolites:
                    self._metabolites.append(m)
                
    @property
    def reactions(self):
        return self._reactions
    
    @property
    def metabolites(self):
        return self._metabolites
    
    @property
    def matrix(self):
        return self._mat
                