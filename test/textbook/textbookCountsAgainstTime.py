from src.generator import FluxGenerator
import matplotlib.pyplot as plt
import cobra.test
from src import models, emailer

# Generate some EFMs
model = cobra.test.create_test_model("textbook")
startReaction = model.reactions.get_by_id("Biomass_Ecoli_core")
include = {startReaction.id: models.FORWARD}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)
countsPath = fg.dumpCountsToFile()
#fg.disableManualStop()
#fg.setMaxCount(1000)
fg.alpha = 2.1
fg.beta = 2.2
fg.gamma = 2
fg.removeDuplicates()

fg.genAll()
 
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
        if len(parts) >= 4:
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
plt.show()
  
#emailer.emailImage("Complete!", body, "Temp.png")


