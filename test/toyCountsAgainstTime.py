from src.generator import FluxGenerator
import matplotlib.pyplot as plt
import cobra
from src import models

# Generate some EFMs
model = cobra.io.load_json_model("toyModel.json")
startReaction = model.reactions.get_by_id("EX_G6P")
include = {startReaction.id: models.FORWARD}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)
fg.setMaxCount(300)
fg.suppressOutput()
#fg.removeDuplicates()

counts = [0]
times = [0]
icounts = [0] # Infeasible count
dcounts = [0] # Duplicate count
ucounts = [0] # Unique count
# Make the extra function retrieve counts
def extra(flux, excludeVal):
    counts.append(len(fg.efmsGenerated))
    ucounts.append(len(fg.uniqueEFMs))
    times.append(fg.getTimeDelta())
    icounts.append(fg.infeasibleCount)
    dcounts.append(fg.duplicateCount)
fg.setExtra(extra)

fg.genAll()

ip, = plt.plot(times, icounts, label="Infeasible", color='red')
dp, = plt.plot(times, dcounts, label="Duplicate", color='green')
up, = plt.plot(times, ucounts, label="Unique", color="blue")
plt.legend(handles=[ip, dp, up])
plt.xlabel("Time")
plt.ylabel("EFM count")
plt.show()