'''
Created on Nov 15, 2016

@author: Lucy
'''

import cobra
import sys
from gurobipy import Model, LinExpr, GRB, tupledict
from incidence import IncidenceMatrix

FORWARD = True
REVERSE = False



def addSteadyStateConstraints(gurobiModel, incidence, externalMetabolites,
                              forwardCoeffs, reverseCoeffs):
    """Adds constraints corresponding to the steady state constraint Cx=0
    
    gurobiModel: The Gurobi model to which the constraints to will be added
    incidence: an IncidenceMatrix object
    externalMetabolites: the metabolites to ignore in steadys tate
    forwardCoeffs: list of Gurobi vars corresponding to the coefficients of the forward reactions
    reverseCoeffs: list of Gurobi vars corresponding to the coefficients of the reverse reactions
    """
    # Cx = 0
    cs = {}
    for m in incidence.metabolites:
        if m.id not in externalMetabolites:
            cs[m.id] = LinExpr()
    for (e, m) in incidence.matrix:
        if m not in externalMetabolites:
            ratio = incidence.matrix[e, m]
            cs[m] += ratio * forwardCoeffs[e]
            if incidence.reactions.get_by_id(e).reversibility:
                cs[m] += -ratio * reverseCoeffs[e]
                
    for c in cs:
        gurobiModel.addConstr(cs[c] == 0)
        
def addReversibilityConstraints(gurobiModel, incidence, forwardCoeffs, reverseCoeffs):
    """Add constraints corresponding to reversibility
    
    gurobiModel: The Gurobi model to which the constraints to will be added
    incidence: an IncidenceMatrix object
    forwardCoeffs: list of Gurobi vars corresponding to the coefficients of the forward reactions
    reverseCoeffs: list of Gurobi vars corresponding to the coefficients of the reverse reactions
    """
    # Allow reversible reactions to have negative coefficients
    # This is done by adding a coefficient for the reverse reaction
    #  and using a SOS to ensure that one is zero
    for r in incidence.reactions:
        if r.reversibility:
            reverseCoeffs[r.id] = gurobiModel.addVar(name=r.id)
            gurobiModel.addSOS(GRB.SOS_TYPE1, [forwardCoeffs[r.id], reverseCoeffs[r.id]])
            
def addIncludedExcludedReactionConstraints(gurobiModel, incidence, reactionsToInclude,
                                           reactionsToExclude, forwardCoeffs, reverseCoeffs):
    """Add constraints corresponding to included and excluded reactions
    
    gurobiModel: The Gurobi model to which the constraints to will be added
    incidence: an IncidenceMatrix object
    reactionsToInclude: dict of reaction IDs and directions, indicating flux should be > 1 or < -1
            in the solution
    reactionsToExclude: list of reaction IDs for which flux = 0 in the solution
    forwardCoeffs: list of Gurobi vars corresponding to the coefficients of the forward reactions
    reverseCoeffs: list of Gurobi vars corresponding to the coefficients of the reverse reactions    
    """
    # Include reactions
    reactions = incidence.reactions
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

def makeGurobiModel(cobraModel, incidence, externalMetabolites, reactionsToInclude, reactionsToExclude):
    gurobiModel = Model(cobraModel.id)
    forwardCoeffs = gurobiModel.addVars([r.id for r in incidence.reactions])
    reverseCoeffs = tupledict()
    
    # ## Constraints
    # Reversibility
    addReversibilityConstraints(gurobiModel, incidence, forwardCoeffs, reverseCoeffs)
    # Cx = 0
    addSteadyStateConstraints(gurobiModel, incidence, externalMetabolites,
                              forwardCoeffs, reverseCoeffs)
    # Include and exclude reactions
    addIncludedExcludedReactionConstraints(gurobiModel, incidence, reactionsToInclude,
                                           reactionsToExclude, forwardCoeffs, reverseCoeffs)
        
    # ## Set objective function
    # Minimise sum of fluxes
    gurobiModel.setObjective(forwardCoeffs.sum() + reverseCoeffs.sum(), GRB.MINIMIZE)
    
    return (gurobiModel, forwardCoeffs, reverseCoeffs)

def extractReactantsAndProducts(flux, reactions):
    reactants = []
    products = []
    
    for f in flux:
        reaction = reactions.get_by_id(f)
        rf = reaction.reactants
        pf = reaction.products
        if flux[f] > 0:
            reactants += [r.id for r in rf]
            products += [p.id for p in pf]
        else:
            reactants += [p.id for p in pf]
            products += [r.id for r in rf]
            
    return (reactants,products)

def main():
    cobraModel = cobra.io.load_json_model("toyModel.json")
    reactionsToInclude = {"Pgi": FORWARD, "Pyk": FORWARD}
    reactionsToExclude = ["Prs_DeoB", "Zwf"]
    externalMetabolites = ["NADP", "NADPH", "CO2", "G6P", "R5Pex", "ATP", "ADP", "Pyr", "NAD", "NADH"]
    
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
    incidence = IncidenceMatrix(reactions, metabolites)
            
    (gurobiModel, forwardCoeffs, reverseCoeffs) = makeGurobiModel(cobraModel, incidence,
                                                              externalMetabolites, reactionsToInclude,
                                                              reactionsToExclude)
            
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
        
        print flux

    (reactants,products) = extractReactantsAndProducts(flux, reactions)

if __name__ == '__main__':
    main()