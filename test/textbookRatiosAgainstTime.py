from src.generator import FluxGenerator
import cobra.test
from src import models
import matplotlib.pyplot as plt

model = cobra.test.create_test_model("textbook")
startReaction = model.reactions.get_by_id("EX_glc__D_e")
include = {startReaction.id: models.REVERSE}
initialExclude = []
fg = FluxGenerator(model, startReaction, include, initialExclude)
fg.suppressOutput()
#fg.useAutoStop()
fg.setMaxCount(400)

times = [0]
iratios = [0] # Infeasible ratios
dratios = [0] # Duplicate ratios
ucounts = [0] # Unique counts

def extra(flux, excludeVal):
    ic = fg.infeasibleCount # Infeasible count
    dc = fg.duplicateCount # Duplicate count
    vc = len(fg.uniqueEFMs) # Valid count
    
    times.append(fg.getTimeDelta())
    iratios.append(float(ic)/vc)
    dratios.append(float(dc)/vc)
    ucounts.append(fg.getUniqueCount())

fg.setExtra(extra)    
fg.genAll()

ip, = plt.plot(times, iratios, label="Infeasible-to-unique", color="red")
dp, = plt.plot(times, dratios, label="Duplicate-to-unique", color="green")
plt.legend(handles=[ip,dp])
plt.xlabel("Time")
plt.ylabel("Ratio")
plt.show()