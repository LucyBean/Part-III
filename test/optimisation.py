from src.generator import FluxGenerator
import cobra
from src import models, emailer
import scipy.optimize
import time

# Load in model
model = cobra.io.read_sbml_model("MODEL1108160000.xml")

# Pick a random reaction to include
#randID = random.randrange(0, len(model.reactions))
#startReaction = model.reactions[randID]
startReaction = model.reactions.get_by_id("ARGSS")

runNum = 0

def runModel(params, *args):
    global runNum
    alpha = params[0]
    beta = params[1]
    print "Run number: {}".format(runNum)
    print "\tAlpha = {:6.4f} Beta = {:6.4f}".format(alpha, beta)
    if alpha <= 0 or beta <= 0:
        return 0
    
    # Try to run the model!
    include = {startReaction.id: models.FORWARD}
    initialExclude = []
    fg = FluxGenerator(model, startReaction, include, initialExclude)
    fg.silent()
    fg.setMaxTime(180)
    fg.suppressOutput()
    fg.removeDuplicates()
    fg.disableManualStop()
    
    def score(r, reacCounts):
        if r not in reacCounts:
            return 1
        fc = reacCounts[r]["f"]
        ic = reacCounts[r]["i"]
        return (fc+1)**alpha /float(fc+ic+1)**beta
    
    fg.score = score
    
    fg.genAll(strategy = 2)
        
    score = len(fg.uniqueEFMs)
    print "\tScore {}".format(score)
    runNum += 1
    return score

startTime = time.time()
res = scipy.optimize.differential_evolution(runModel, [(0.1,2.0), (0.1,2.0)], maxiter=4, popsize=10, disp=True)

deltaTime = time.time() - startTime
textResults = """Results:
{res!s}

Time taken: {time}""".format(res = res, time = deltaTime, disp=True)

print textResults
emailer.emailText("Optimisation complete", textResults)

