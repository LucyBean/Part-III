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
startReaction = model.reactions.get_by_id("FACOAL161t2pp")


print "Including reaction", startReaction.id

# Try to run the model!
include = {startReaction.id: models.FORWARD}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)
#fg.useAutoStop(ratio=2)
#fg.setMaxTime(3600)
fg.disableManualStop()
fg.setMaxCount(500)
fg.removeDuplicates()
countsPath = fg.dumpCountsToFile()
fg.alpha = 6.4
fg.beta = 6.6
fg.gamma = 10

fg.genAll()

if len(fg.efmsGenerated) > 1:
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
            # parts[1] shows total number generated (dc+uc)
            icounts.append(parts[2])
            dcounts.append(parts[3])
            ucounts.append(parts[4])
        
    plt.clf()
    ip, = plt.plot(times, icounts, label="Infeasible", color='red')
    dp, = plt.plot(times, dcounts, label="Duplicate", color='green')
    up, = plt.plot(times, ucounts, label="Unique", color="blue")
    plt.legend(handles=[ip, dp, up])
    plt.xlabel("Time")
    plt.ylabel("EFM count")
    plt.savefig("Temp.png")
    emailer.emailImage("Complete!", body, "Temp.png")
