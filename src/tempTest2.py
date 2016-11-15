'''
Created on Nov 15, 2016

@author: Lucy
'''
obj = LinExpr()
for m in metabolites:
    if m.id not in externalMetabolites:
        cs[m.id] = LinExpr()
for (e, m) in incidence:
    ratio = incidence[e, m]
    if m in externalMetabolites:
        if ratio > 0:
            print forwardCoeffs[e]
            obj += ratio * forwardCoeffs[e]
        elif reactions.get_by_id(e).reversibility:
            print reverseCoeffs[e]
            obj += ratio * reverseCoeffs[e]