from src.generator import FluxGenerator
import matplotlib.pyplot as plt
import cobra
from src import models
import time
import random

# Load in model
model = cobra.io.read_sbml_model("MODEL1108160000.xml")

# Pick a random reaction to include
#randID = random.randrange(0, len(model.reactions))
#startReaction = model.reactions[randID]
startReaction = model.reactions.get_by_id("ARGSS")

print "Including reaction", startReaction.id

# Try to run the model!
include = {startReaction.id: models.FORWARD}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)
#fg.useAutoStop(ratio=2)
fg.setMaxTime(300)
#fg.setMaxCount(100)
fg.suppressOutput()
fg.removeDuplicates()
 
counts = [0]
times = [0]
icounts = [0] # Infeasible count
dcounts = [0] # Duplicate count
ucounts = [0] # Unique count
mcounts = [0] # Minimal count
 
startTime = time.time()
# Make the extra function check the count
def extra(flux, excludeVal, **kwargs):
    counts.append(len(fg.efmsGenerated))
    ucounts.append(len(fg.uniqueEFMs))
    times.append(fg.getTimeDelta())
    icounts.append(fg.infeasibleCount)
    dcounts.append(fg.duplicateCount)
    mcounts.append(len(fg.getMinimalEFMs()))
fg.setExtra(extra)
 
fg.genAll()
 
if len(fg.efmsGenerated) > 0:
    # Show the plot if some were generated
    ip, = plt.plot(times, icounts, label="Infeasible", color='red')
    dp, = plt.plot(times, dcounts, label="Duplicate", color='green')
    up, = plt.plot(times, ucounts, label="Unique", color="blue")
    mp, = plt.plot(times, mcounts, label="Minimal", color="purple")
    plt.legend(handles=[ip, dp, up, mp])
    plt.xlabel("Time")
    plt.ylabel("EFM count")
    plt.show()
