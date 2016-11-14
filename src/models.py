'''
Created on Nov 12, 2016

@author: Lucy
'''
from gurobipy import Model, LinExpr, GRB, tupledict
import escher
import sys

FORWARD = True
REVERSE = False

def process(model,
                 reactionsToInclude=[],
                 reactionsToExclude=[],
                 ignoredMetabolites=[]):
    """Process a model, producing a dictionary containing all included reactions
    and their fluxes.
    
    model: A COBRA model to be used.
    reactionsToExclude: A list of reactions whose fluxes should be zero.
    reactionsToInclude: A dict containing reaction/direction pairs. Direction should be
        set to either models.FORWARD or models.BACKWARD
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
    forwardCoeffs = model.addVars([r.id for r in reactions])
    reverseCoeffs = tupledict()
    
    
    # ## Constraints
    # Allow reversible reactions to have negative coefficients
    # This is done by adding a coefficient for the reverse reaction
    #  and using a SOS to ensure that one is zero
    for r in reactions:
        if r.reversibility:
            reverseCoeffs[r.id] = model.addVar(name=r.id)
            model.addSOS(GRB.SOS_TYPE1, [forwardCoeffs[r.id], reverseCoeffs[r.id]])
    # Cx = 0
    cs = {}
    for m in metabolites:
        if m.id not in ignoredMetabolites:
            cs[m.id] = LinExpr()
    for (e, m) in incidence:
        if m not in ignoredMetabolites:
            ratio = incidence[e, m]
            cs[m] += ratio * forwardCoeffs[e]
            if reactions.get_by_id(e).reversibility:
                cs[m] += -ratio * reverseCoeffs[e]
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
        
    ### Set objective function
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
                    raise BaseException("Forward and reverse for " + c  + " used in solution.")
                flux[c] = -v
        
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
