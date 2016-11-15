'''
Created on Nov 15, 2016

@author: Lucy
'''

import cobra
import sys
from gurobipy import Model, LinExpr, GRB, tupledict

FORWARD = True
REVERSE = False

if __name__ == '__main__':
    ### Calculating p-semiflows
    gurobiModel = cobra.io.load_json_model("toyModel.json")
    metabolitesToInclude=["Ru5P"]
    
    metabolites = gurobiModel.metabolites
    reactions = gurobiModel.reactions
    incidence = {}
    for r in reactions:
        for m in r.metabolites:
            incidence[r.id, m.id] = r.get_coefficient(m)
            
    # ## Make gurobi model
    gurobiModel = Model(gurobiModel.id)
    coeffs = gurobiModel.addVars([m.id for m in metabolites])
    
    
    # ## Constraints
    # yC = 0
    cs = {}
    for r in reactions:
        cs[r.id] = LinExpr()
    for (e, m) in incidence:
        ratio = incidence[e, m]
        cs[e] += ratio * coeffs[m]
    for c in cs:
        gurobiModel.addConstr(cs[c] == 0)
        
    # TODO: Include metabolites
    for m in metabolitesToInclude:
        # Check if the metabolite name is valid
        if m in [a.id for a in metabolites]:
            gurobiModel.addConstr(coeffs[m] >= 1)
        else:
            sys.stderr.write("Attempting to include unknown metabolite " + m)
        
    # ## Set objective function
    # Minimise sum of fluxes
    gurobiModel.setObjective(coeffs.sum(), GRB.MINIMIZE)
    
    gurobiModel.optimize()
    
    if gurobiModel.status == GRB.OPTIMAL:
        flux = {}
        for c in coeffs:
            v = coeffs[c].x
            if v > 0.01:
                flux[c] = v
        
        print flux