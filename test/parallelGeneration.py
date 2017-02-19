from mpi4py import MPI
import models
import cobra.test
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
batch_size = 1

def genResult(exclude):
	"""Produces a flux from an exclude set"""
	flux = models.findEFM(model, include, exclude, 0)
	return flux

if rank == 0:
	# Master
	
	model = cobra.test.create_test_model("textbook")
	startReaction = "EX_glc__D_e"
	include = {startReaction: models.REVERSE}
	initialExclude = set([])

	alpha = 6.4
	beta = 6.6
	gamma = 10

	infeasibleCount = 0
	duplicateCount = 0

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
		global reacCounts
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
		global toTryNum
		
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
	
	def processResults(resultBatch):
		global infeasibleCount
		global duplicateCount
		
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
			
	def makeBatch():
		global toTryNum
		global batch_size
		global reacScores
		
		excludeSets = []
		while toTryNum > 0 and len(excludeSets) < batch_size:
			next = pickNext(toTry, reacScores)
			if next is None:
				break
			else:
				nextReaction, exclude = next
			excludeSets.append(next)
		return excludeSets
		
	startTime = time.time()
	
	# Init queues
	flux = genResult(initialExclude)
	tried.append(initialExclude)
	results.append(flux)
	processResults([(flux, startReaction, initialExclude)])
	
	assert toTryNum > 0
	
	suspendedProcs = []
	workingProcs = size
	
	# Send a batch to each slave
	for i in range(1,size):
		excludeSets = makeBatch()
			
		if len(excludeSets) > 0:
			# Send the next batch
			data = {"type":"vals", "vals":excludeSets}
			comm.send(data, dest=i)
		else:
			# Make this process wait
			suspendedProcs.append(i)
			workingProcs -= 1
		
	# Get results from each slave
	while workingProcs > 1:
		data = comm.recv(source = MPI.ANY_SOURCE)
		slaveRank = data["from"]
		# Get result and record slave as 'waiting'
		if data["type"] == "res":
			res = data["res"]
			processResults(res)
			suspendedProcs.append(slaveRank)
			workingProcs -= 1
			
		# Try to unsuspend processes
		unsuspended = []
		for i in suspendedProcs:
			if checkStop():
				print "Stopped"
				break
				
			excludeSets = makeBatch()				
			if len(excludeSets) > 0:
				# Send the next batch
				data = {"type":"vals", "vals":excludeSets}
				comm.send(data, dest=i)
				unsuspended.append(i)
				workingProcs += 1
			
		suspendedProcs = [p for p in suspendedProcs if p not in unsuspended]
			
	# Stop all suspended procs
	for i in suspendedProcs:
		comm.send({"type":"stop"},i)
	
	deltaTime = time.time() - startTime
	print """#Results: {}
Time taken: {}""".format(len(results), deltaTime)
	
else:
	# Slave
	
	# Receive instruction
	while True:
		data = comm.recv(source=0)
		
		if data["type"] == "vals":
			# Check if values were received
			# Get values received
			excludeSets = data["vals"]
			
			# Compute results
			results = []
			for (nr,es) in excludeSets:
				flux = genResult(es)
				results.append((flux, nr, es))
			
			# Send results
			data = {"from":rank, "type":"res", "res":results}
			comm.send(data, dest=0)
			
		if data["type"] == "stop":
			# Check if told to stop
			break
			
print rank, "has finished"
		