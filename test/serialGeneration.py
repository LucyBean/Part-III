from src import models
import cobra.test
import time

model = cobra.test.create_test_model("textbook")
startReaction = "EX_glc__D_e"
include = {startReaction: models.REVERSE}
initialExclude = set([])

alpha = 6.4
beta = 6.6
gamma = 10

infeasibleCount = 0
duplicateCount = 0

solvingTime = 0
processingTime = 0
pickingTime = 0

maxAttempts = 1000
maxTime = -1

results = []
resultsByNames = []
toTry = {}
toTryNum = 0
tried = []

reacCounts = {}
reacScores = {}

stopping = False

def addCount(name, count):
	""" Adds a count to this reaction
		count = {"f":1} will add 1 to feasible count
	"""
	if name not in reacCounts:
		reacCounts[name] = {}
		reacCounts[name]["d"] = 0
		reacCounts[name]["u"] = 0
		reacCounts[name]["i"] = 0
	for k in count:
		reacCounts[name][k] += count[k]

def score(name, counts):
	if name not in counts:
		return 100
	uc = counts[name]["u"]
	dc = counts[name]["d"]
	ic = counts[name]["i"]
	
	return float(uc+1)**alpha / (dc*gamma + ic + 1) ** beta

def fluxToNames(flux):
	return set(flux.keys())
	
def addToTryList(reaction, exclude):
	global toTryNum
	if reaction not in toTry:
		toTry[reaction] = []
	toTry[reaction].append(exclude)
	toTryNum += 1
	
def pickNext(toTry, reacScores):
	"""Given a list of (nextReaction, exclude) types, returns the next tuple to try.
	Splices them out of the toTry list"""
	global pickingTime
	global toTryNum
	
	pickingStart = time.time()
	def getScore(nextReaction):
		if nextReaction in reacScores:
			return reacScores[nextReaction]
		else:
			return 20
	
	# Sort all the next possible reactions by score
	nextReacs = sorted(toTry.keys(), key = lambda k: getScore(k))
	nextReacs = [nr for nr in nextReacs if len(toTry[nr]) > 0]
	i = 0
	# Find the next exclude set for this reaction that has not been tried
	while toTryNum > 0:
		while len(toTry[nextReacs[i]]) == 0:
			# Find first non empty
			i += 1
			if i == len(toTry):
				break
		if i == len(toTry):
			break
		excludeSets = toTry[nextReacs[i]]
		es = excludeSets[0]
		excludeSets = toTry[nextReacs[i]] = excludeSets[1:]
		toTryNum -= 1
		if es not in tried:
			pickingTime += time.time() - pickingStart
			return nextReacs[i], es
			
def output(string):
	print string
	
def getTotalSolves():
	return len(results) + duplicateCount + infeasibleCount
		
def stop(reason):
	global stopping
	if reason is not None:
		output("Stopping: {}".format(reason))
	stopping = True
		
def checkStop():
	"""Indicates whether the generation loop should stop"""
	if toTryNum == 0:
		stop("Explored all EFMs")
	if (maxAttempts > 0	and getTotalSolves() > maxAttempts):
		stop("Tried maximum number of EFMs")
	if (maxTime > 0	and deltaTime > maxTime):
		stop("Reached maximum execution time")
		
def getResults():
	s = """
Generated:
{unique:>6} unique EFMs
{dupl:>6} duplicate EFMs
{infs:>6} infeasible EFMs""".format(unique = len(results),
							 dupl = duplicateCount,
							 infs = infeasibleCount)
	return s
		
def printResults(self):
	self.output(self.getResults())

def genNextExcludeSets(exclude, flux):
	"""Produces all other exclude sets to test from this flux"""
	fluxNames = fluxToNames(flux)
	knockOutAble = [r for r in fluxNames if r not in include and r not in exclude]
	for r in knockOutAble:
		nexclude = exclude | set([r])
		addToTryList(r, nexclude)

def genResult(exclude):
	"""Produces a flux from an exclude set"""
	global solvingTime
	solveStart = time.time()
	flux = models.findEFM(model, include, exclude, 0)
	solvingTime += time.time() - solveStart
	return flux

def processResults(resultBatch):
	global infeasibleCount
	global duplicateCount
	global processingTime
	
	processStart = time.time()
	for (flux, nextReaction, exclude) in resultBatch:
		if flux is None:
			# No solution found
			addCount(nextReaction, {"i":1})
			infeasibleCount += 1
			continue
			
		names = fluxToNames(flux)
		if names in resultsByNames:
			# Duplicate result
			addCount(nextReaction, {"d":1})
			duplicateCount += 1
			continue
			
		else:
			# Unique result
			addCount(nextReaction, {"u":1})
			results.append(flux)
			resultsByNames.append(names)
			genNextExcludeSets(exclude, flux)
			
		# Update the score
		reacScores[nextReaction] = score(nextReaction, reacCounts)
		processingTime += time.time() - processStart
		
startTime = time.time()
deltaTime = 0
exclude = initialExclude
nextReaction = startReaction
while True:
	# Generate a result
	tried.append(exclude)
	flux = genResult(exclude)
	processResults([(flux, nextReaction, exclude)])
	
	# Pick the next to try
	n = pickNext(toTry, reacScores)
	if n is None:
		stop("Explored all EFMs")
	else:
		nextReaction, exclude = n
	
	deltaTime = time.time() - startTime
	
	print " "*70, "\r",
	checkStop()
	if stopping:
		break
	else:
		print "Time: {:6.2f} Count: {:7} Infeasible: {:7} Total solves: {:7}\r".format(deltaTime, len(results), infeasibleCount, getTotalSolves()),

output(getResults())
output("Time taken: {:.4f}".format(deltaTime))
		