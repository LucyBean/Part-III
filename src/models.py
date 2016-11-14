'''
Created on Nov 12, 2016

@author: Lucy
'''
from gurobipy import Model, LinExpr, GRB
import escher

def process(model,
                 reactionsToInclude=[],
                 reactionsToExclude=[],
                 ignoredMetabolites=[]):
    """Process a model, producing a dictionary containing all included reactions
    and their fluxes.
    
    model: A COBRA model to be used.
    reactionsToExclude: A list of reactions whose fluxes should be zero.
    reactionsToInclude: A list of reactions whose fluxes should be >= 1.
    ignoredMetabolites: Metabolites which are ignored when finding a solution.    
    """
    
    # ## Extract the needed data from the model
    metabolites = model.metabolites
    reactions = model.reactions
    incidence = {}
    for r in reactions:
        for m in r.metabolites:
            incidence[r.id, m.id] = r.get_coefficient(m)
            
    # ## Make gurobi model
    model = Model("textbook")
    coeffs = model.addVars([r.id for r in reactions])
    
    
    # ## Constraints
    # Allow reversible reactions to have negative coefficients
    for r in reactions:
        if r.reversibility:
            coeffs[r.id].lb = -GRB.INFINITY
    # Cx = 0
    cs = {}
    for m in metabolites:
        if m.id not in ignoredMetabolites:
            cs[m.id] = LinExpr()
    for (e, m) in incidence:
        if m not in ignoredMetabolites:
            ratio = incidence[e, m]
            cs[m] += ratio * coeffs[e]
    for c in cs:
        model.addConstr(cs[c] == 0)
    # Include and exclude reactions
    for r in reactionsToInclude:
        model.addConstr(coeffs[r] >= 1)
    for r in reactionsToExclude:
        model.addConstr(coeffs[r] == 0)
        
    ### Set objective function
    model.setObjective(coeffs.sum(), GRB.MINIMIZE)
            
    model.optimize()
    
    if model.status == GRB.OPTIMAL:
        flux = {}
        for c in coeffs:
            v = coeffs[c].x
            if v != 0:
                flux[c] = v
        
        return flux
    

def display(map_name=None, map_json=None, reaction_data=[]):
    b = escher.Builder(map_name=map_name, map_json=map_json,
                       reaction_data=reaction_data,
                       # color and size according to the absolute value
                   reaction_styles=['color', 'size', 'abs', 'text'],
                   # change the default colors
                   reaction_scale=[{'type': 'min', 'color': '#00cc00', 'size': 4},
                                   {'type': 'mean', 'color': '#0000dd', 'size': 20},
                                   {'type': 'max', 'color': '#ff0000', 'size': 40}])
    b.display_in_browser(scroll_behavior="zoom")
