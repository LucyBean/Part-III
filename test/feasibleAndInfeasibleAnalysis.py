import json
import matplotlib.pyplot as plt

# Find scores for each reaction
def score(fc, ic):
    threshold = 10
    if (fc + ic) < threshold:
        return 1
    else:
        return ((fc+1) / float(fc+ic+1))**0.4 * (fc+1)**0.5

if "data" not in locals():
    with open("ARGSS feasible and infeasible.json", "r") as f:
        s = f.read()
        data = json.loads(s)
    
fdata = data["feasibles"]
idata = data["infeasibles"]

reacCounts = {}

# For each reaction, get number of times
#   it makes feasible and infeasible exclude set
for f in fdata:
    reaction = f[0]
    if reaction not in reacCounts:
        reacCounts[reaction] = {}
        reacCounts[reaction]["f"] = 0
        reacCounts[reaction]["i"] = 0
    reacCounts[reaction]["f"] += 1
    
for i in idata:
    reaction = i[0]
    if reaction not in reacCounts:
        reacCounts[reaction] = {}
        reacCounts[reaction]["f"] = 0
        reacCounts[reaction]["i"] = 0
    reacCounts[reaction]["i"] += 1
    
# Find totals and ratio
for r in reacCounts:
    reacCounts[r]["total"] = reacCounts[r]["f"] + reacCounts[r]["i"]
    reacCounts[r]["ratio"] = reacCounts[r]["i"] / float(reacCounts[r]["total"])
    reacCounts[r]["score"] = score(reacCounts[r]["f"], reacCounts[r]["i"])
    
# Find distribution of lengths of exclusion sets
fLens = [len(f[1]) for f in fdata]
iLens = [len(i[1]) for i in idata]
fx = [i/float(len(fLens)) for i in range(len(fLens))]
ix = [i/float(len(iLens)) for i in range(len(iLens))]

# Find "maybe bads"
maybeBads = [r for r in reacCounts if reacCounts[r]["ratio"] > 0 and reacCounts[r]["ratio"] < 1]
ratios = [reacCounts[r]["ratio"] for r in reacCounts]
totals = [reacCounts[r]["total"] for r in reacCounts]
    
def plotRatios():
    # plot the distribution of ratios
    plt.plot(sorted(ratios))
    plt.ylabel("Infeasible ratio")
    plt.xlabel("#Reaction")
    plt.show()

def plotHist():
    plt.hist(ratios, 20)
    plt.ylabel("Density")
    plt.xlabel("Infeasible ratio")
    plt.show()
    
def plotTotals():
    plt.plot(sorted(totals))
    plt.ylabel("Total number of EFMs with this reaction")
    plt.xlabel("#Reaction")
    plt.show()
    
def plotTotals2():
    goodTotals = [reacCounts[r]["i"] for r in reacCounts if reacCounts[r]["ratio"] < 0.1]
    badTotals = [reacCounts[r]["i"] for r in reacCounts if reacCounts[r]["ratio"] > 0.9]
    g, = plt.plot(sorted(goodTotals), label="Ratio < 0.1")
    b, = plt.plot(sorted(badTotals), label="Ratio > 0.9")
    plt.ylabel("Total number of EFMs with this reaction")
    plt.xlabel("#Reaction")
    plt.legend(handles=[g,b])
    plt.show()
    
def plotLens():
    f, = plt.plot(fx, sorted(fLens), label="Feasible")
    i, = plt.plot(ix, sorted(iLens), label="Infeasible")
    plt.ylabel("Length of EFM")
    plt.xlabel("%EFM")
    plt.legend(handles=[f,i])
    plt.show()
    
def plotScores():
    rs = [reacCounts[r]["ratio"] for r in reacCounts]
    ss = [reacCounts[r]["score"] for r in reacCounts]
    plt.plot(rs, ss, 'o')
    plt.xlabel("Ratio")
    plt.ylabel("Score")
    plt.show()
    
    
    