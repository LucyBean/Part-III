'''
Created on Nov 12, 2016

@author: Lucy
'''
from gurobipy import Model, LinExpr, GRB, tupledict
from incidence import IncidenceMatrix
from src import helpers
import sys
import json

FORWARD = True
REVERSE = False
logToConsole = True

def addSteadyStateConstraints(gurobiModel, incidence, forwardCoeffs, reverseCoeffs):
    """Adds constraints corresponding to the steady state constraint Cx=0
    
    gurobiModel: The Gurobi model to which the constraints to will be added
    incidence: an IncidenceMatrix object
    forwardCoeffs: list of Gurobi vars corresponding to the coefficients of the forward reactions
    reverseCoeffs: list of Gurobi vars corresponding to the coefficients of the reverse reactions
    """
    # Cx = 0
    cs = {}
    for m in incidence.metabolites:
        cs[m.id] = LinExpr()
    for (e, m) in incidence.matrix:
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

def makeGurobiModel(cobraModel, reactionsToInclude, reactionsToExclude):
    """Makes a Gurobi model from the given Cobra model, ready to be optimized.
    The objective function used will minimise the the sum of absolute flux
    through the reactions.    
    """
    # ## Extract the needed data from the cobra model
    incidence = IncidenceMatrix(cobraModel.reactions)
    
    gurobiModel = Model(cobraModel.id)
    forwardCoeffs = gurobiModel.addVars([r.id for r in incidence.reactions])
    reverseCoeffs = tupledict()
    
    # ## Constraints
    # Reversibility
    addReversibilityConstraints(gurobiModel, incidence, forwardCoeffs, reverseCoeffs)
    # Cx = 0
    addSteadyStateConstraints(gurobiModel, incidence, forwardCoeffs, reverseCoeffs)
    # Include and exclude reactions
    addIncludedExcludedReactionConstraints(gurobiModel, incidence, reactionsToInclude,
                                           reactionsToExclude, forwardCoeffs, reverseCoeffs)
    # Minimum flux sum constraint
    gurobiModel.addConstr(forwardCoeffs.sum() + reverseCoeffs.sum() >= 1)
        
    # ## Set objective function
    # Minimise sum of fluxes
    gurobiModel.setObjective(forwardCoeffs.sum() + reverseCoeffs.sum(), GRB.MINIMIZE)
    return (gurobiModel, forwardCoeffs, reverseCoeffs)

def findTerminalReactantsAndProducts(flux, reactions):
    """Find the metabolites that are at the start and end of this pathway.
    
    flux: A dict of fluxes, as given by models.findEFM
    reactions: A list of the reactions, given by cobraModel.reactions
    
    return: (terminalReactants, terminalProducts)
    """
    externalReactants = []
    externalProducts = []
    
    # Find all external reactions which have a flux through them
    externalReactions = [r.id for r in reactions if r.id[0:3] == "EX_"]
    for ep in externalReactions:
        if ep in flux:
            # If the flux through this reaction is positive then the metabolite
            #   is being consumed by the reaction
            f = flux[ep]
            
            # Find the metabolite involved in this reaction
            epr = reactions.get_by_id(ep)
            if epr.products != []:
                metabolite = epr.products[0]
            else:
                metabolite = epr.reactants[0]
                
            # Append to reactants/products
            if f > 0:
                externalReactants.append(metabolite)
            else:
                externalProducts.append(metabolite)
    
    terminalProducts = []
    terminalReactants = []
    
    # Find the final product
    for ep in externalProducts:
        # find all the non-external reactions that are used in solution
        myReactions = [r for r in ep.reactions if r.id[0:3] != "EX_" and r.id in flux]
        allExternalReactions = True
    
        for mr in myReactions:
            if flux[mr.id] > 0:
                products = mr.products
            else:
                products = mr.reactants
                
            # Check if all the products are external
            allExternalMetabolites = True
            for p in products:
                if p not in externalProducts:
                    allExternalMetabolites = False
                    break
                
            if not allExternalMetabolites:
                allExternalReactions = False
            
        if allExternalReactions:    
            terminalProducts.append(ep)
                
    # Find the initial reactant
    for ep in externalReactants:
        # find all the non-external reactions that are used in solution
        myReactions = [r for r in ep.reactions if r.id[0:3] != "EX_" and r.id in flux]
        allExternalReactions = True
        
        for mr in myReactions:
            if flux[mr.id] > 0:
                reactants = mr.reactants
            else:
                reactants = mr.products
                
            # Check if all the products are external
            allExternalMetabolites = True
            for p in reactants:
                if p not in externalReactants:
                    allExternalMetabolites = False
                    break
                
            if not allExternalMetabolites:
                allExternalReactions = False
                
        if allExternalReactions:
            terminalReactants.append(ep)
                
    return (terminalReactants, terminalProducts)

def findEFM(cobraModel,
                 reactionsToInclude=[],
                 reactionsToExclude=[],
                 logToConsole = 1):
    """Process a gurobiModel, producing a dictionary containing all included reactions
    and their fluxes. Uses an objective function that minimises the sum of the
    fluxes.
    
    gurobiModel: A COBRA gurobiModel to be used.
    reactionsToExclude: A list of reactions whose fluxes should be zero.
    reactionsToInclude: A dict containing reaction/direction pairs. Direction should be
        set to either models.FORWARD or models.BACKWARD  
    """
            
    (gurobiModel, forwardCoeffs, reverseCoeffs) = makeGurobiModel(cobraModel, reactionsToInclude,
                                                              reactionsToExclude)
    gurobiModel.params.logtoConsole = logToConsole
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
                    print forwardCoeffs
                    print reverseCoeffs
                    raise BaseException("Forward and reverse for " + c + " used in solution.")
                flux[c] = -v
        
        return flux
    
def findProducts(cobraModel, startID):
    metabolites = cobraModel.metabolites
    reactions = cobraModel.reactions
    if startID not in metabolites:
        sys.stderr.write("Attempting to start from unknown metabolite " + startID)
        return
    startingMetabolite = metabolites.get_by_id(startID)
    # This will get all possible starting reactions. This excludes all reactions
    #  that are marked as external (start with "EX_")
    possibleStartReactions = [r for r in list(startingMetabolite.reactions) if r.id[0:3] != "EX_"]
    products = {}
    
    for startReaction in possibleStartReactions:
        reactionsToInclude = {startReaction.id: FORWARD}
        reactionsToExclude = []
        
        # Discover up to 10 reactions per initial start reaction.
        for _ in range(10):
            flux = findEFM(cobraModel, reactionsToInclude, reactionsToExclude)            
            
            # If no solutions remain then break out of loop
            if flux is None:
                break
            
            (_, externalProducts) = findTerminalReactantsAndProducts(flux, reactions)
            # Remove the starting metabolite if it appears in this list
            externalProducts = [p for p in externalProducts if p.id != startID]
        
            # Find the reaction that immediately leads to production of terminal product
            if externalProducts != []:
                terminalProduct = externalProducts[0]
                rs = list(terminalProduct.reactions)
                # Take the first reaction in this list
                knockOut = [r for r in rs if r.id in flux and r.id[0:3] != "EX_"][0]
            # The EMF produced is a loop
            else:
                terminalProduct = startingMetabolite
                rs = list(terminalProduct.reactions)
                # Find all the reactions that the starting metabolite is involved in
                mightKnockOut = [r for r in rs if r.id in flux and r.id[0:3] != "EX_"]
                # Find all of the ones that produce the starting metabolite
                knockOutList = []
                for r in mightKnockOut:
                    if flux[r.id] > 0 and startingMetabolite in r.products:
                        knockOutList.append(r)
                    elif flux[r.id] < 0 and startingMetabolite in r.reactants:
                        knockOutList.append(r)
                knockOut = knockOutList[0]
                
            # Display output
            if logToConsole:
                print "Found flux to", terminalProduct.id
                print "Include:", reactionsToInclude
                print "Exclude:", reactionsToExclude
                helpers.printFlux(flux)
                print "\n"
            
            # Knock out reaction that produces the terminal product
            reactionsToExclude.append(knockOut.id)
            
            # Add this possible product/pathway to the products dict
            tpid = terminalProduct.id
            if tpid not in products:
                products[tpid] = {}
                products[tpid]["desc"] = json.dumps("Include: " + str(reactionsToInclude) + "<br/>Exclude: " + str(reactionsToInclude))
                products[tpid]["fluxes"] = []
            products[tpid]["fluxes"].append(flux)
            
    return products
    
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


