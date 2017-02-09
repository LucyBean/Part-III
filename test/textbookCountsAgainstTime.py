from src.generator import FluxGenerator
import matplotlib.pyplot as plt
import cobra.test
from src import models, emailer
import time

# Generate some EFMs
model = cobra.test.create_test_model("textbook")
startReaction = model.reactions.get_by_id("EX_glc__D_e")
include = {startReaction.id: models.REVERSE}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)
fg.setMaxTime(60)
fg.setMaxCount(1000)
fg.suppressOutput()
#fg.useManualInput()
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

fg.genAll(2)

body = fg.getConfig() + "\n\n---Results---\n" + fg.getResults()
#emailer.emailText("Complete!","The run has completed.")

plt.clf()
ip, = plt.plot(times, icounts, label="Infeasible", color='red')
dp, = plt.plot(times, dcounts, label="Duplicate", color='green')
up, = plt.plot(times, ucounts, label="Unique", color="blue")
mp, = plt.plot(times, mcounts, label="Minimal", color="purple")
plt.legend(handles=[ip, dp, up, mp])
plt.xlabel("Time")
plt.ylabel("EFM count")
plt.savefig("Temp.png")

emailer.emailImage("Complete!", body, "Temp.png")


