'''
Created on Nov 12, 2016

@author: Lucy
'''
import cobra.test
import json
from gurobipy import Model, LinExpr, GRB
textbook_model = cobra.test.create_test_model("textbook")


### Extract model data
# Metabolites and reaction names
metabolites = textbook_model.metabolites
reactions = textbook_model.reactions
# Incidence matrix
incidence = {}
for r in textbook_model.reactions:
    for m in r.metabolites:
        incidence[r.id,m.id] = r.get_coefficient(m)
        
reactionsToInclude = ["SUCDi"]
reactionsToExclude = ["FRD7"]
ignoredMetabolites = ["adp_c", "atp_c", "coa_c", "h2o_c", "h_c", "h_e",
                      "nad_c", "nadh_c", "nadp_c", "nadph_c", "pep_c"]
    
    
### Make gurobi model
model = Model("textbook")
coeffs = model.addVars([r.id for r in reactions])


### Constraints
# Allow reversible reactions to have negative coefficients
for r in reactions:
    if r.reversibility:
        coeffs[r.id].lb = -GRB.INFINITY
        
# Cx = 0
cs = {}
for m in metabolites:
    if m.id not in ignoredMetabolites:
        cs[m.id] = LinExpr()
for (e,m) in incidence:
    if m not in ignoredMetabolites:
        ratio = incidence[e,m]
        cs[m] += ratio * coeffs[e]
for c in cs:
    model.addConstr(cs[c] == 0)
    
# Include and exclude reactions
for r in reactionsToInclude:
    model.addConstr(coeffs[r] >= 1)
for r in reactionsToExclude:
    model.addConstr(coeffs[r] == 0)
        
model.optimize()

if model.status == GRB.OPTIMAL:
    flux = {}
    for c in coeffs:
        v = coeffs[c].x
        if v != 0:
            flux[c] = v
            print c,v
    
    with open("reactions.json", "w") as f:
        json.dump(flux, f)