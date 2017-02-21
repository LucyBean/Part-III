from src.generator import FluxGenerator
import matplotlib.pyplot as plt
import cobra
from src import models, emailer
#import random

# Load in model
model = cobra.io.read_sbml_model("MODEL1108160000.xml")

# Pick a random reaction to include
#randID = random.randrange(0, len(model.reactions))
#startReaction = model.reactions[randID]
#startReaction = model.reactions.get_by_id("ARGSS")
#startReaction = model.reactions.get_by_id("FACOAL161t2pp")
startReaction = model.reactions.get_by_id("Ec_biomass_iJO1366_WT_53p95M")


print "Including reaction", startReaction.id

# Try to run the model!
include = {startReaction.id: models.FORWARD}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)
#fg.useAutoStop(ratio=2)
#fg.setMaxTime(3600)
#fg.verboseOutput()
#fg.useManualInput()
#fg.disableManualStop()
#fg.setMaxCount(500)
fg.removeDuplicates()
countsPath = fg.dumpCountsToFile()
fg.alpha = 2.1
fg.beta = 2.2
fg.gamma = 2

fg.genAll()

if len(fg.uniqueEFMs) > 1:
    body = fg.getConfig() + "\n\n---Results---\n" + fg.getResults()
    
    times = []
    icounts = []
    dcounts = []
    ucounts = []
    # Load the counts from the file
    with open(countsPath, "r") as countsFile:
        # Skip header
        countsFile.readline()
        # Parse files
        for line in countsFile:
            parts = line.split(",")
            times.append(parts[0])
            icounts.append(parts[1])
            dcounts.append(parts[2])
            ucounts.append(parts[3])
        
    plt.clf()
    ip, = plt.plot(times, icounts, label="Infeasible", color='red')
    dp, = plt.plot(times, dcounts, label="Duplicate", color='green')
    up, = plt.plot(times, ucounts, label="Unique", color="blue")
    plt.legend(handles=[ip, dp, up])
    plt.xlabel("Time")
    plt.ylabel("EFM count")
    plt.savefig("Temp.png")
    emailer.emailImage("Complete!", body, "Temp.png")
