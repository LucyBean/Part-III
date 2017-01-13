from src.generator import FluxGenerator
import cobra
from src import models

model = cobra.io.load_json_model("toyModel.json")
startReaction = model.reactions.get_by_id("EX_G6P")
include = {startReaction.id: models.FORWARD, "EX_Pyr": models.REVERSE}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)

print "Use manual input? (y/n)",
if (raw_input() == "y"):
    fg.useManualInput()
    
print "Max number to generate: (enter -1 for no limit)",
maxCount = int(raw_input())
if maxCount > 0:
    fg.setMaxCount(maxCount)
    
fg.genAll()