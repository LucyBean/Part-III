from src.generator import FluxGenerator
import cobra
from src import models, helpers, display

# Generate some EFMs
model = cobra.io.load_json_model("toyModel.json")
startReaction = model.reactions.get_by_id("EX_G6P")
include = {startReaction.id: models.FORWARD}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)
fg.setMaxCount(100)
fg.suppressOutput()
fg.genAll()
minimalEFMs = fg.getMinimalEFMs()
    
print "Minimal EFMs:"
for m in minimalEFMs:
    print "\t", m

print "Input 'r' to show results",
showResults = (raw_input() == 'r')

if showResults:
    # Display minimal
    toDisplay = {}
    for m in minimalEFMs:
        reactionIDs = helpers.reactionValToNames(m, fg.reactionVals)
        flux = {}
        for r in reactionIDs:
            flux[r] = 1
        toDisplay[m] = {"fluxes":[flux]}
        
    display.displayAll(map_json = "toyModelMap.json", toDisplay=toDisplay, title="Minimal EFMs")
