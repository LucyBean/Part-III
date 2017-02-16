from src.generator import FluxGenerator
import cobra
from src import models

model = cobra.io.load_json_model("toyModel.json")
startReaction = model.reactions.get_by_id("EX_G6P")
include = {startReaction.id: models.FORWARD}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)
fg.dumpCountsToFile()
fg.disableManualStop()
fg.alpha = 2.6
fg.beta = 2.7
fg.genAll(1)