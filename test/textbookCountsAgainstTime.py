from src.generator import FluxGenerator
import matplotlib.pyplot as plt
import cobra.test
from src import models
import time

# Generate some EFMs
model = cobra.test.create_test_model("textbook")
startReaction = model.reactions.get_by_id("EX_glc__D_e")
include = {startReaction.id: models.REVERSE}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)
fg.setMaxTime(150)
fg.setMaxCount(1000)
fg.suppressOutput()
#fg.removeDuplicates()

counts = [0]
times = [0]
icounts = [0] # Infeasible count
dcounts = [0] # Duplicate count
ucounts = [0] # Unique count
mcounts = [0] # Minimal count

startTime = time.time()
# Make the extra function check the count
def extra(flux, excludeVal):
    count = len(fg.efmsGenerated)
    prevCount = counts[-1]
    if (count != prevCount):
        counts.append(count)
        ucounts.append(len(set(fg.efmsGenerated)))
        times.append(fg.getTimeDelta())
        icounts.append(fg.infeasibleCount)
        dcounts.append(fg.duplicateCount)
        mcounts.append(len(fg.getMinimalEFMs()))
fg.setExtra(extra)

fg.genAll()

#plt.plot(times, icounts, 'ro', times, dcounts, 'go', times, ucounts, 'bo')
ip, = plt.plot(times, icounts, label="Infeasible", color='red')
dp, = plt.plot(times, dcounts, label="Duplicate", color='green')
up, = plt.plot(times, ucounts, label="Unique", color="blue")
mp, = plt.plot(times, mcounts, label="Minimal", color="purple")
plt.legend(handles=[ip, dp, up, mp])
plt.xlabel("Time")
plt.ylabel("EFM count")
plt.show()