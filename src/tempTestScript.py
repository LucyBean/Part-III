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
    model = cobra.io.load_json_model("toyModel.json")
    reactionsToInclude={"Pgi":FORWARD}
    reactionsToExclude=[]
    externalMetabolites=["NADP", "NADPH", "CO2", "G6P", "R5Pex", "ATP", "ADP", "Pyr", "NAD", "NADH"]
    
    metabolites = model.metabolites
    reactions = model.reactions
    incidence = {}
    for r in reactions:
        for m in r.metabolites:
            incidence[r.id, m.id] = r.get_coefficient(m)
            
    # ## Make gurobi model
    model = Model("textbook")
    forwardCoeffs = model.addVars([r.id for r in reactions])
    reverseCoeffs = tupledict()
    
    
    # ## Constraints
    # Allow reversible reactions to have negative coefficients
    # This is done by adding a coefficient for the reverse reaction
    #  and using a SOS to ensure that one is zero
    for r in reactions:
        if r.reversibility:
            reverseCoeffs[r.id] = model.addVar()
            model.addSOS(GRB.SOS_TYPE1, [forwardCoeffs[r.id], reverseCoeffs[r.id]])
    # Cx = 0
    cs = {}
    obj = LinExpr()
    for m in metabolites:
        if m.id not in externalMetabolites:
            cs[m.id] = LinExpr()
    for (e, m) in incidence:
        ratio = incidence[e, m]
        if m not in externalMetabolites:
            cs[m] += ratio * forwardCoeffs[e]
            if reactions.get_by_id(e).reversibility:
                cs[m] += -ratio * reverseCoeffs[e]
        else:
            if ratio > 0:
                obj += ratio * forwardCoeffs[e]
            elif reactions.get_by_id(e).reversibility:
                obj += -ratio * reverseCoeffs[e]
    for c in cs:
        model.addConstr(cs[c] == 0)
    # Include reactions
    for r in reactionsToInclude:
        # Check if the reaction name is valid
        if r in [a.id for a in reactions]:
            direction = reactionsToInclude[r]
            if direction == FORWARD:
                model.addConstr(forwardCoeffs[r] >= 1)
            # If adding a reverse, check the reaction is reversible
            elif reactions.get_by_id(r).reversibility:
                model.addConstr(reverseCoeffs[r] >= 1)
            else:
                sys.stderr.write("Attempting to force reverse for irreversible reaction " + r + "\n")
        else:
            sys.stderr.write("Attempting to include unknown reaction:" + r + "\n")
    # Exclude reactions
    for r in reactionsToExclude:
        if r in [a.id for a in reactions]:
            model.addConstr(forwardCoeffs[r] == 0)
            if reactions.get_by_id(r).reversibility:
                model.addConstr(reverseCoeffs[r] == 0)
        else:
            sys.stderr.write("Attempting to exclude unknown reaction:" + r + "\n")
        
    # ## Set objective function
    # Minimise sum of fluxes
    model.setObjective(forwardCoeffs.sum() + reverseCoeffs.sum(), GRB.MINIMIZE)
    
    model.optimize()
    
    if model.status == GRB.OPTIMAL:
        flux = {}
        for c in forwardCoeffs:
            v = forwardCoeffs[c].x
            if v > 0.01:
                flux[c] = v
        for c in reverseCoeffs:
            v = reverseCoeffs[c].x
            if v > 0.01:
                if c in flux.keys():
                    raise BaseException("Forward and reverse for " + c + " used in solution.")
                flux[c] = -v
        
        print flux