'''
Created on Nov 11, 2016

@author: Lucy
'''
from gurobipy import Model, GRB, LinExpr

metabolites = ["Ru5P","FP2","F6P","GAP","R5P"]
revReactions = ["pgi","fba","rpi","rpe"]
irrevReactions = ["gap","zwf","pfk","fbp","prs"]
reactions = revReactions + irrevReactions
reactionsToInclude = ["prs","gap"]
reactionsToExclude = ["zwf"]

#Incidence matrix as a dict
incidence = {
    ("pgi","F6P"):1,
    ("fba","FP2"):-1, ("fba","GAP"):2,
    ("rpi","Ru5P"):-1, ("rpi","R5P"):1,
    ("rpe","Ru5P"):-2, ("rpe","F6P"):2, ("rpe","GAP"):1, ("rpe","R5P"):-1,
    ("gap","GAP"):-1,
    ("zwf","Ru5P"):1,
    ("pfk","FP2"):1, ("pfk","F6P"):-1,
    ("fbp","FP2"):-1, ("fbp","F6P"):1,
    ("prs","R5P"):-1
}


#Make models
model = Model("incidence")
coeffs = model.addVars(reactions)


### Constraints
# Allow reversible reactions to have negative co-efficients
for r in revReactions:
    coeffs[r].lb = -GRB.INFINITY
    # This does not work, see updated "models.py" in src
    
# Cx = 0
cs = {}
for m in metabolites:
    cs[m] = LinExpr()
for (e,m) in incidence:
    ratio = incidence[e,m]
    cs[m] += ratio * coeffs[e]
for c in cs:
    model.addConstr(cs[c] == 0)
    
# Include and exclude reactions
for r in reactionsToInclude:
    model.addConstr(coeffs[r] >= 1)
for r in reactionsToExclude:
    model.addConstr(coeffs[r] == 0)


### Optimize
model.optimize()
for c in coeffs:
    if coeffs[c].x > 0:
        print c , coeffs[c].x