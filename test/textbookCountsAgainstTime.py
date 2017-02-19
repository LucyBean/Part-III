from src.generator import FluxGenerator
import matplotlib.pyplot as plt
import cobra.test
from src import models, emailer

# Generate some EFMs
model = cobra.test.create_test_model("textbook")
startReaction = model.reactions.get_by_id("EX_glc__D_e")
include = {startReaction.id: models.REVERSE}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)
countsPath = fg.dumpCountsToFile()
fg.disableManualStop()
fg.setMaxCount(581)
fg.alpha = 6.4
fg.beta = 6.6
fg.gamma = 10
fg.removeDuplicates()

fg.genAll(1)
 
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
        if len(parts) >= 5:
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
plt.show()
  
#emailer.emailImage("Complete!", body, "Temp.png")


