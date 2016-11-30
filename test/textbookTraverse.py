'''
Created on Nov 30, 2016

@author: Lucy
'''
import cobra.test
from src import models

if __name__ == '__main__':
    if "flux" not in locals():
        model = cobra.test.create_test_model("textbook")
        include = {"EX_glc__D_e":models.REVERSE}
        exclude = []
        flux = models.findEFM(model, include, exclude)
        
    reactions = model.reactions
    startReaction = reactions.get_by_id("EX_glc__D_e")
    reactionsToExplore = [startReaction]
    exploredReactions = []
    
    while reactionsToExplore != []:
        # Remove the reaction at the tail of the list
        thisReaction = reactionsToExplore[-1]
        reactionsToExplore = reactionsToExplore[:-1]
        exploredReactions.append(thisReaction)
        print "Exploring", thisReaction.id
        
        # Find the products of this reaction
        if flux[thisReaction.id] > 0:
            products = thisReaction.products
        else:
            products = thisReaction.reactants
        
        # Find the next metabolite in the pathway
        if len(products) == 0:
            # This is an external reaction so no more reactions available
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
        nextReactions = [r for r in nextMetabolite.reactions if r.id in flux and r not in exploredReactions]
        # Ensure it is a reactant in this reaction
        for r in nextReactions:
            if flux[r.id] > 0 and nextMetabolite in r.products:
                nextReactions.remove(r)
            elif flux[r.id] < 0 and nextMetabolite in r.reactants:
                nextReactions.remove(r)
                
        if len(nextReactions) > 1:
            print "\tbranch"
        # Append these reactions to reactionsToExplore
        reactionsToExplore = reactionsToExplore + nextReactions
        