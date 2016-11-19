'''
Created on Nov 12, 2016

@author: Lucy
'''
from gurobipy import Model, LinExpr, GRB, tupledict
import escher
import sys

FORWARD = True
REVERSE = False


def findEFM(cobraModel,
                 reactionsToInclude=[],
                 reactionsToExclude=[],
                 externalMetabolites=[]):
    """Process a gurobiModel, producing a dictionary containing all included reactions
    and their fluxes. Uses an objective function that minimises the sum of the
    fluxes.
    
    gurobiModel: A COBRA gurobiModel to be used.
    reactionsToExclude: A list of reactions whose fluxes should be zero.
    reactionsToInclude: A dict containing reaction/direction pairs. Direction should be
        set to either models.FORWARD or models.BACKWARD
    externalMetabolites: Metabolites which are ignored when finding a solution.    
    """
    
    # ## Extract the needed data from the cobra model
    metabolites = cobraModel.metabolites
    reactions = cobraModel.reactions
    incidence = {}
    for r in reactions:
        for m in r.metabolites:
            incidence[r.id, m.id] = r.get_coefficient(m)
            
    # ## Make gurobi model
    gurobiModel = Model(cobraModel.id)
    forwardCoeffs = gurobiModel.addVars([r.id for r in reactions])
    reverseCoeffs = tupledict()
    
    
    # ## Constraints
    # Allow reversible reactions to have negative coefficients
    # This is done by adding a coefficient for the reverse reaction
    #  and using a SOS to ensure that one is zero
    for r in reactions:
        if r.reversibility:
            reverseCoeffs[r.id] = gurobiModel.addVar(name=r.id)
            gurobiModel.addSOS(GRB.SOS_TYPE1, [forwardCoeffs[r.id], reverseCoeffs[r.id]])
    # Cx = 0
    cs = {}
    for m in metabolites:
        if m.id not in externalMetabolites:
            cs[m.id] = LinExpr()
    for (e, m) in incidence:
        if m not in externalMetabolites:
            ratio = incidence[e, m]
            cs[m] += ratio * forwardCoeffs[e]
            if reactions.get_by_id(e).reversibility:
                cs[m] += -ratio * reverseCoeffs[e]
    for c in cs:
        gurobiModel.addConstr(cs[c] == 0)
    # Include reactions
    for r in reactionsToInclude:
        # Check if the reaction name is valid
        if r in [a.id for a in reactions]:
            direction = reactionsToInclude[r]
            if direction == FORWARD:
                gurobiModel.addConstr(forwardCoeffs[r] >= 1)
            # If adding a reverse, check the reaction is reversible
            elif reactions.get_by_id(r).reversibility:
                gurobiModel.addConstr(reverseCoeffs[r] >= 1)
            else:
                sys.stderr.write("Attempting to force reverse for irreversible reaction " + r + "\n")
        else:
            sys.stderr.write("Attempting to include unknown reaction:" + r + "\n")
    # Exclude reactions
    for r in reactionsToExclude:
        if r in [a.id for a in reactions]:
            gurobiModel.addConstr(forwardCoeffs[r] == 0)
            if reactions.get_by_id(r).reversibility:
                gurobiModel.addConstr(reverseCoeffs[r] == 0)
        else:
            sys.stderr.write("Attempting to exclude unknown reaction:" + r + "\n")
        
    # ## Set objective function
    # Minimise sum of fluxes
    gurobiModel.setObjective(forwardCoeffs.sum() + reverseCoeffs.sum(), GRB.MINIMIZE)
            
    gurobiModel.optimize()
    
    
    if gurobiModel.status == GRB.OPTIMAL:
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
        
        return flux
      
def findEFM_alt_obj(cobraModel,
                 reactionsToInclude=[],
                 reactionsToExclude=[],
                 externalMetabolites=[]):
    """Process a gurobiModel, producing a dictionary containing all included reactions
    and their fluxes. Uses an objective function that minimises the sum of the
    fluxes.
    
    gurobiModel: A COBRA gurobiModel to be used.
    reactionsToExclude: A list of reactions whose fluxes should be zero.
    reactionsToInclude: A dict containing reaction/direction pairs. Direction should be
        set to either models.FORWARD or models.BACKWARD
    externalMetabolites: Metabolites which are ignored when finding a solution.    
    """
    
    # ## Extract the needed data from the cobra model
    metabolites = cobraModel.metabolites
    reactions = cobraModel.reactions
    incidence = {}
    for r in reactions:
        for m in r.metabolites:
            incidence[r.id, m.id] = r.get_coefficient(m)
            
    # ## Make gurobi model
    gurobiModel = Model(cobraModel.id)
    forwardCoeffs = gurobiModel.addVars([r.id for r in reactions])
    reverseCoeffs = tupledict()
    
    
    # ## Constraints
    # Allow reversible reactions to have negative coefficients
    # This is done by adding a coefficient for the reverse reaction
    #  and using a SOS to ensure that one is zero
    for r in reactions:
        if r.reversibility:
            reverseCoeffs[r.id] = gurobiModel.addVar()
            gurobiModel.addSOS(GRB.SOS_TYPE1, [forwardCoeffs[r.id], reverseCoeffs[r.id]])
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
        gurobiModel.addConstr(cs[c] == 0)
    # Include reactions
    for r in reactionsToInclude:
        # Check if the reaction name is valid
        if r in [a.id for a in reactions]:
            direction = reactionsToInclude[r]
            if direction == FORWARD:
                gurobiModel.addConstr(forwardCoeffs[r] >= 1)
            # If adding a reverse, check the reaction is reversible
            elif reactions.get_by_id(r).reversibility:
                gurobiModel.addConstr(reverseCoeffs[r] >= 1)
            else:
                sys.stderr.write("Attempting to force reverse for irreversible reaction " + r + "\n")
        else:
            sys.stderr.write("Attempting to include unknown reaction:" + r + "\n")
    # Exclude reactions
    for r in reactionsToExclude:
        if r in [a.id for a in reactions]:
            gurobiModel.addConstr(forwardCoeffs[r] == 0)
            if reactions.get_by_id(r).reversibility:
                gurobiModel.addConstr(reverseCoeffs[r] == 0)
        else:
            sys.stderr.write("Attempting to exclude unknown reaction:" + r + "\n")
        
    # ## Set objective function
    # Minimise sum of fluxes
    gurobiModel.setObjective(forwardCoeffs.sum() + reverseCoeffs.sum(), GRB.MINIMIZE)
    
    gurobiModel.optimize()
    
    if gurobiModel.status == GRB.OPTIMAL:
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
        
        return flux
    
def findPsemi(cobraModel, metabolitesToInclude):
    metabolites = cobraModel.metabolites
    reactions = cobraModel.reactions
    incidence = {}
    for r in reactions:
        for m in r.metabolites:
            incidence[r.id, m.id] = r.get_coefficient(m)
            
    # ## Make gurobi model
    gurobiModel = Model(cobraModel.id)
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
        
        return flux

def displayEFM(map_name=None, map_json=None, reaction_data=[]):
    b = escher.Builder(map_name=map_name, map_json=map_json,
                       reaction_data=reaction_data,
                       # color and size according to the absolute value
                   reaction_styles=['color', 'size', 'abs', 'text'],
                   # change the default colors
                   reaction_scale=[{'type': 'min', 'color': '#00cc00', 'size': 4},
                                   {'type': 'mean', 'color': '#0000dd', 'size': 20},
                                   {'type': 'max', 'color': '#ff0000', 'size': 40}])
    b.display_in_browser(scroll_behavior="zoom")
    
def displayPsemi(map_name=None, map_json=None, metabolite_data=[]):
    b = escher.Builder(map_name=map_name, map_json=map_json,
                       metabolite_data=metabolite_data,
                       # color and size according to the absolute value
                   reaction_styles=['color', 'size', 'abs', 'text'],
                   # change the default colors
                   reaction_scale=[{'type': 'min', 'color': '#00cc00', 'size': 4},
                                   {'type': 'mean', 'color': '#0000dd', 'size': 20},
                                   {'type': 'max', 'color': '#ff0000', 'size': 40}])
    b.display_in_browser(scroll_behavior="zoom")
