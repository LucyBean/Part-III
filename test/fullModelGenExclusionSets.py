from src.generator import FluxGenerator
import matplotlib.pyplot as plt
import cobra
from src import models
import time
import random

# Load in model
model = cobra.io.read_sbml_model("MODEL1108160000.xml")

# Pick a random reaction to include
randID = random.randrange(0, len(model.reactions))
startReaction = model.reactions[randID]
#startReaction = model.reactions.get_by_id("ARGSS")

print "Including reaction", startReaction.id

# Try to run the model!
include = {startReaction.id: models.FORWARD}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)
#fg.useAutoStop(ratio=2)
fg.disableManualStop()
fg.setMaxTime(600)
#fg.setMaxCount(100)
fg.suppressOutput()
#fg.removeDuplicates()
 
feasibles = []
infeasibles = []

def extra(flux, excludeSet, **kwargs):
    feasibles.append((kwargs["newest"], excludeSet))
fg.setExtra(extra)

def infeasibleExtra(newest, excludeSet, **kwargs):
    infeasibles.append((newest, excludeSet))
fg.setInfeasibleExtra(infeasibleExtra)
 
fg.genAll()