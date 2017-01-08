import cobra
from src import models, display
import json
import time

def reactionNamesToVal(reactionIDs):
    val = 0
    for r in reactionIDs:
        val += reactionVals[r]
    return val

model = cobra.io.load_json_model("toyModel.json")
reactions = model.reactions
startReaction = reactions.get_by_id("EX_G6P")
include = {startReaction.id: models.FORWARD, "EX_Pyr": models.REVERSE}
exclude = []
flux = models.findEFM(model, include, exclude, 0)

reactionVals = {}
nextVal = 1
for r in reactions:
    reactionVals[r.id] = nextVal
    nextVal *= 2

efmsToExplore = [(flux, exclude)]
efmsGenerated = [reactionNamesToVal(flux.keys())]
testedExclusions = [reactionNamesToVal(exclude)]
displayOutput = {}
displayOutput[0] = {}
displayOutput[0]["desc"] = json.dumps("Include: " + str(include) + "<br/>Exclude: " + str(exclude))
displayOutput[0]["fluxes"] = [flux]
# This will keep track of which exclusion sets produced each EFM
efmsGeneratedByExclusionSet = {}

count = 0
infeasibleCount = 0
duplicateCount = 0

print "Use manual input? (y/n)",
useManual = (raw_input() == "y")

print "Max number to generate: (enter -1 for no limit)",
maxCount = int(raw_input())

startTime = time.time()

inp = ""

while inp != "q":
    if efmsToExplore == []:
        print "\nExplored all EFMs"
        break
    if maxCount > 0 and len(efmsGenerated) > maxCount:
        print "\nGenerated max number of EFMs"
        break
    print "Exploring EFM", count
    efm = efmsToExplore[-1]
    efmsToExplore = efmsToExplore[:-1]
    flux = efm[0]
    exclude = efm[1]
    count += 1
    
    knockOutAble = [r for r in reactions if r.id not in include and r.id not in exclude]
    for r in knockOutAble:
        # Create a new exclusion set with this reaction and the exclusion set used to generate this
        nexclude = exclude + [r.id]
        nexclude.sort()
        nexcludeVal = reactionNamesToVal(nexclude)
        print "Next exclusion set", nexcludeVal
        
        if useManual:
            inp = raw_input() # Pause to allow quitting or looking at input
            if inp == "q":
                break
        
        if nexcludeVal in testedExclusions: # See if this exclusion set has been tested
            print "\tDuplicate exclusion set."
            continue # Skip this exclusion set
        
        testedExclusions.append(nexcludeVal)
        # Generate the EFM
        nflux = models.findEFM(model, include, nexclude, 0)
        
        if nflux is not None:
            # Find the value for this efm
            efmReactions = nflux.keys()
            efmReactionsVal = reactionNamesToVal(efmReactions)
            print "\tFlux generated", efmReactionsVal
            
            efmsToExplore.append((nflux, nexclude))
            efmsGenerated.append(efmReactionsVal)
            # Check whether max number of EFMs have been generated
            if maxCount > 0 and len(efmsGenerated) >= maxCount:
                break
            # See if this has been generated before and increment count if so
            if efmReactionsVal not in efmsGeneratedByExclusionSet:
                efmsGeneratedByExclusionSet[efmReactionsVal] = []
                # Save the display
                i = efmReactionsVal
                displayOutput[i] = {}
                displayOutput[i]["desc"] = json.dumps("Include: " + str(include) + "<br/>Exclude: " + str(nexclude))
                displayOutput[i]["fluxes"] = [nflux]
            else:
                duplicateCount += 1
                print "\tEFM has been generated before"
            efmsGeneratedByExclusionSet[efmReactionsVal].append(efmReactionsVal)
            
        else:
            infeasibleCount += 1
            print "\tInfeasible"
            
# Optionally display results
print "Enter r to display results",
inp = raw_input()

endTime = time.time()
print "Generated", len(efmsGenerated), "in", (endTime - startTime), "s"
print "Generated", infeasibleCount, "infeasible EFMs"
print "Generated", duplicateCount, "duplicate EFMs"
    
if inp == "r":
    display.displayAll(map_json="toyModelMap.json", toDisplay=displayOutput, title="All knockouts")
            