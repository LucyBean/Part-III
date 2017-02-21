'''
Created on Nov 30, 2016

@author: Lucy
'''
import cobra.test
from src import models

if "flux" not in locals():
    model = cobra.test.create_test_model("textbook")
    include = {"EX_glc__D_e":models.REVERSE}
    exclude = []
    flux = models.findEFM(model, include, exclude, 0)
        
reactions = model.reactions
startReaction = reactions.get_by_id("EX_glc__D_e")

reactionsToExplore = [startReaction]
exploredReactions = []
attemptedKnockOuts = []
finalProducts = []

efmsToExplore = []

while reactionsToExplore != []:
    # Remove the reaction at the tail of the list
    thisReaction = reactionsToExplore[-1]
    reactionsToExplore = reactionsToExplore[:-1]
    if thisReaction in exploredReactions:
        continue
    exploredReactions.append(thisReaction)
    print "Exploring", thisReaction.id
    
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
            nexclude = exclude + [nr.id]
            nflux = models.findEFM(model, include, nexclude, 0)
            
            print "\t", nexclude
            if nflux is not None:
                print "\tFlux generated"
                efmsToExplore.append(nflux)
            else:
                print "\tEFM is none."
            
            print ""
        