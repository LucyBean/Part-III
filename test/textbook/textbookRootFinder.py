from src.generator import FluxGenerator
import cobra.test
from src import models, helpers
import random

model = cobra.test.create_test_model("textbook")

interestingIDs = ["PGL","THD2","ATPM"]

#===============================================================================
# numTrials = 10
# for _ in range(numTrials):
#     startReactionID = random.randrange(0, len(model.reactions))
#     startReaction = model.reactions[startReactionID]
#===============================================================================
for startReactionID in interestingIDs:
    startReaction = model.reactions.get_by_id(startReactionID)
    print "Start reaction", startReaction.id
    
    include = {startReaction.id: models.FORWARD}
    initialExclude = []
    fg = FluxGenerator(model, startReaction, include, initialExclude)
    reactionVals = helpers.genReactionVals(model.reactions)
    
    # This will keep track of which exclusion sets produced each EFM
    efmsGeneratedByExclusionSet = {}
    roots = {}
    # Results table
    results = {"c":0, "n":0}
    def extra(flux, excludeSet):
        global roots
        global efmsGeneratedByExclusionSet
        global results
        # Get s
        efmReactionsVal = helpers.reactionNamesToVal(flux.keys(), reactionVals)
        
        # Initialise dict if neccesary
        if efmReactionsVal not in efmsGeneratedByExclusionSet:
            print "Unique:", excludeSet
            efmsGeneratedByExclusionSet[efmReactionsVal] = []
        
        # Check whether roots need to be initialised
        
        if efmReactionsVal not in roots:
            nonEmptyESs = [ne for ne in efmsGeneratedByExclusionSet[efmReactionsVal] if len(ne) > 0]
            if len(nonEmptyESs) > 0:
                thisLen = len(excludeSet)
                minLen = len(min(nonEmptyESs))
                if thisLen > minLen:
                    roots[efmReactionsVal] = set()
                    for es in efmsGeneratedByExclusionSet[efmReactionsVal]:
                        for k in es:
                            roots[efmReactionsVal] |= set([k])
                    print "\tRoots set for", efmReactionsVal
                
        # If the roots have been initialised
        else:
            # Check if this contains a root
            myRoots = roots[efmReactionsVal]
            containsRoot = False
            for root in myRoots:
                if root in excludeSet:
                    containsRoot = True
                    break
            # Record in table
            k = "n"
            if containsRoot:
                k = "c"
            results[k] += 1
                
        # Store the generating exclusion set for this flux
        efmsGeneratedByExclusionSet[efmReactionsVal].append(excludeSet)
        
    fg.setExtra(extra)
    fg.setMaxCount(1000)
    #fg.useManualInput()
    fg.suppressOutput()
    fg.genAll()
    
    print "\t", results