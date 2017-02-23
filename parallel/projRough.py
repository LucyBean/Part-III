from mpi4py import MPI
import random
from time import sleep, time
import numpy as np
import cobra
import models
from generator import FluxGenerator
import datetime
import sys
import json

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

batch_size = 10
max_results = 200

results = []

# Clear file
myFilePath = "out/{}.txt".format(rank)

sys.stderr = open("out/{}err.txt".format(rank), "w")

with open(myFilePath,"w") as f:
	f.write("Rank: {}\n".format(rank))
	f.write("Time: {}\n".format(datetime.datetime.now()))
	
# Set up flux generator
# Load in model
model = cobra.io.read_sbml_model("MODEL1108160000.xml")
startReaction = model.reactions.get_by_id("Ec_biomass_iJO1366_WT_53p95M")
include = {startReaction.id: models.FORWARD}
initialExclude = set([])
fg = FluxGenerator(model, startReaction, include, initialExclude)
fg.silent()
fg.removeDuplicates()
fg.alpha = 2.1
fg.beta = 2.2
fg.gamma = 2

def output(string, end="\n"):
	with open(myFilePath,"a",0) as f:
		f.write(string + end)

def genNextExcludeSets(exclude, flux):
	return fg.genNextExcludeSets(exclude, flux)
			
def findSolution(excludeSet):
	return fg.findSolution(excludeSet)
	
def processResult(flux, nextReaction, exclude):
	return fg.processResult(flux, nextReaction, exclude)

def processResults(batch):
	for (flux, nextReaction, exclude) in batch:
		if flux is not None:
			results.append(flux)
		processResult(flux, nextReaction, exclude)
		
def makeBatch():
	# Makes a batch of values to try and splices them out
	batch = []
	for _ in range(batch_size):
		n = pickNext()
		if n is None:
			break
		batch.append(n)
	return batch
	
def pickNext():
	return fg.pickNext()
		
def send(data, **kwargs):
	assert "dest" in kwargs
	ds = str(data)
	if len(ds) > 53:
		ds = ds[:50] + "..."
	output("Sent {} to {}".format(ds, kwargs["dest"]))
	
def recv(**kwargs):
	data = comm.recv(**kwargs)
	assert "source" in kwargs
	ds = str(data)
	if len(ds) > 53:
		ds = ds[:50] + "..."
	output("Received {} from {}".format(ds, kwargs["source"]))
	return data
				
if rank == 0:
	startTime = time()
	# Init queues
	initResult = findSolution(initialExclude)
	results.append(initResult)
	
	genNextExcludeSets(initialExclude, initResult)
	
	suspendedProcs = []
	workingProcs = size
	
	stop = False
	
	output("Initialised")
	
	# Send a batch to each slave
	for i in range(1,size):
		vals = makeBatch()
		if len(vals) > 0:
			data = {"type":"vals", "vals":vals}
			send(data, dest=i)
		else:
			suspendedProcs.append(i)
			output("{} suspended".format(i))
			workingProcs -= 1
		
	# Get results from each slave
	while not stop:
		data = recv(source = MPI.ANY_SOURCE)
		slaveRank = data["from"]
		# Get result and record slave as waiting
		if data["type"] == "res":
			res = data["res"]
			processResults(res)
			suspendedProcs.append(slaveRank)
			workingProcs -= 1
			
		# Try to unsuspend processes
		unsuspended = []
		for i in suspendedProcs:
			vals = makeBatch()
			if len(vals) == 0:
				break
			data = {"type":"vals", "vals":vals}
			send(data, dest=i)
			unsuspended.append(i)
			workingProcs += 1
		suspendedProcs = [p for p in suspendedProcs if p not in unsuspended]
		
		if workingProcs == 1:
			break
		if len(results) > max_results:
			break
		
	output("Stopping all suspended processes {}".format(suspendedProcs))
	# Stop all suspended procs
	for i in suspendedProcs:
		send({"type":"stop"}, dest=i)
	output("Working procs: {}, Suspended procs: {}".format(workingProcs, suspendedProcs))
	
	deltaTime = time() - startTime
	
	def prettyString(flux):
		s = "{\n"
		for k in flux:
			s += "\t{}: {:.5f}\n".format(k, flux[k])
		s += "}"
		return s
		
	with open("out/prettyResults.txt","w") as f:
		f.write("Results")
		for r in results:
			f.write(prettyString(r))
			
	with open("out/results.json","w") as f:
		s = json.dumps(results)
		f.write(s)
			
	output("""	
ToTry: {}

Time taken: {}

#Results: {}""".format(fg.toTryNum, deltaTime, len(results)))
	print "Time taken: {}\nNum results: {}".format(deltaTime, len(results))
		
else:
	# Receive instruction
	while True:
		data = recv(source=0)
		
		if data["type"] == "vals":
			# Check if values were received
			# Get values received
			vals = data["vals"]
			output("Received {}".format(vals))
			
			# Compute results
			results = []
			for (nextReaction, excludeSet) in vals:
				r = findSolution(excludeSet)
				results.append((r, nextReaction, excludeSet))
			
			# Send results
			data = {"from":rank, "type":"res", "res":results}
			send(data, dest=0)
			
		if data["type"] == "stop":
			# Check if told to stop
			output("Stopped")
			break
		

		