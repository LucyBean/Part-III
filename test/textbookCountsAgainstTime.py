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
fg.setMaxCount(500)
fg.alpha = 6.4
fg.beta = 6.6
fg.gamma = 10
fg.removeDuplicates()


def pickNext(efmsToExplore, reacScores):
    """Given a list of (nextReaction, exclude) tuples, picks the next reaction
    
    Reassign this to change the search strategy. By default, picks a random reaction weighted
    by score.
    
    Note: score was calculated when this tuple was added to the list and may not necessarily
    be the same as that reaction's current score
    
    Return: index of the next tuple to use"""
    
    def getScore(nextReaction):
        if nextReaction in reacScores:
            return reacScores[nextReaction]
        else:
            return 20
    # Find index of max score
    nextReacs = [efm[0] for efm in efmsToExplore]
    max_i = max(range(len(nextReacs)), key = lambda k: getScore(nextReacs[k]))
    return max_i

def score(nextReaction, reacCounts):
    """Calculates the score of a reaction given the reacCounts
    
    reacCounts should be a dictionary which contains unique, feasible, and infeasible counts
    for each reaction."""
    uc = reacCounts[nextReaction]["u"]
    fc = reacCounts[nextReaction]["f"]
    dc = fc - uc # Duplicate count
    ic = reacCounts[nextReaction]["i"]
    
    score = float(uc+1)**fg.alpha / (dc*fg.gamma +ic+1)**fg.beta
    return score

fg.score = score
fg.pickNext = pickNext
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


