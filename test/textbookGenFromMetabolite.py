import cobra.test
from src import models, display
import json
import time

startTime = time.time()

model = cobra.test.create_test_model("textbook")
include = {"EX_glc__D_e":models.REVERSE}
exclude = []
exclude.sort()
flux = models.findEFM(model, include, exclude, 0)
        
reactions = model.reactions
startReaction = reactions.get_by_id("EX_glc__D_e")

displayOutput = {}
displayOutput[0] = {}
displayOutput[0]["desc"] = json.dumps("Include: " + str(include) + "<br/>Exclude: " + str(exclude))
displayOutput[0]["fluxes"] = [flux]

efmsToExplore = [(flux, exclude)]
efmsGenerated = [flux]
efmsGenerated = [flux.keys()]
testedExclusions = [exclude]
count = 0
infeasibleCount = 0
duplicateCount = 0
infeasibleTime = 0
duplicateTime = 0

while efmsToExplore != []:
    print "Exploring efm", count
    count += 1
    # Remove the efm at the tail of the list
    efm = efmsToExplore[-1]
    flux = efm[0]
    efmsToExplore = efmsToExplore[:-1]
    
    reactionsToExplore = [startReaction]
    exploredReactions = []
    attemptedKnockOuts = []
    finalProducts = []
    
    while reactionsToExplore != []:
        # Remove the reaction at the tail of the list
        thisReaction = reactionsToExplore[-1]
        reactionsToExplore = reactionsToExplore[:-1]
        if thisReaction in exploredReactions:
            continue
        exploredReactions.append(thisReaction)
        
        # Find the products of this reaction
        if flux[thisReaction.id] > 0:
            reactants = thisReaction.reactants
            products = thisReaction.products
        else:
            reactants = thisReaction.products
            products = thisReaction.reactants
        
        # Find the next metabolite in the pathway
        if len(products) == 0:
            # This is an external reaction so no more reactions available
            finalProducts += reactants
            continue
        elif len(products) == 1:
            # There is only one product of this reaction
            nextMetabolite = products[0]
        else:
            # If there are multiple metabolites to choose from
            #   choose the one that is involved in the fewest
            #   reactions
            nextMetabolite = products[0]
            for p in products[1:]:
                if len(p.reactions) < len(nextMetabolite.reactions):
                    nextMetabolite = p
            
        # Find the reactions in which nextMetabolite is involved
        nextReactions = [r for r in nextMetabolite.reactions if r.id in flux]
        # Ensure it is a reactant in this reaction
        for r in nextReactions:
            if flux[r.id] > 0 and nextMetabolite in r.products:
                nextReactions.remove(r)
            elif flux[r.id] < 0 and nextMetabolite in r.reactants:
                nextReactions.remove(r)
                
        # Ensure that nextReactions does not contain any reactions that will be considered
        # Sort to explore the reaction with fewest metabolites next
        nextReactions = [r for r in nextReactions if r not in reactionsToExplore]
        nextReactions = sorted(nextReactions, key= lambda r: len(r.metabolites), reverse=True)
                
        # Append these reactions to reactionsToExplore
        reactionsToExplore = reactionsToExplore + nextReactions
        
        # If this is a branch then need to try to generate EFMs with parts of the
        #   branch knocked out
        if len(nextReactions) > 1:
            # Try to knock out each reaction in turn
            for nr in nextReactions:
                if nr in attemptedKnockOuts:
                    continue
                attemptedKnockOuts.append(nr)
                nexclude = efm[1] + [nr.id]
                nexclude.sort()
                if nexclude in testedExclusions:
                    print "\tDuplicate exclusion set."
                    continue
                testedExclusions.append(nexclude)
                t1 = time.time()
                nflux = models.findEFM(model, include, nexclude, 0)
                
                if nflux is not None and nflux != {}:
                    efmReactions = sorted(nflux.keys())
                    if efmReactions not in efmsGenerated:
                        print "\tNew flux generated:", nexclude
                        efmsToExplore.append((nflux, nexclude))
                        efmsGenerated.append(nflux)
                        efmsGenerated.append(efmReactions)
                        i = len(displayOutput)
                        displayOutput[i] = {}
                        displayOutput[i]["desc"] = json.dumps("Include: " + str(include) + "<br/>Exclude: " + str(nexclude))
                        displayOutput[i]["fluxes"] = [nflux]
                    else:
                        print "\tDuplicate flux generated", nexclude
                        duplicateCount += 1
                        duplicateTime += time.time() - t1
                else:
                    print "\tInfeasible", nexclude
                    infeasibleCount += 1
                    infeasibleTime += time.time() - t1

elapsedTime = time.time() - startTime
print "Generated", len(efmsGenerated), "EFMs."
print infeasibleCount, "infeasible EFMs tried"
print duplicateCount, "duplicate EFMs generated"
print "in", elapsedTime, "seconds"
print "Spent", infeasibleTime, "on trying infeasible EFMs"
print "Spent", duplicateTime, "on generating duplicate EFMs"

display.displayAll(map_json="e_coli_core.json", toDisplay=displayOutput, title="EFMs from glc__D_e")
