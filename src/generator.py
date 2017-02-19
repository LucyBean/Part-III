'''
Created on Jan 12, 2017

@author: Lucy
'''

from __future__ import print_function
from src import models
import time
import threading
from msvcrt import getch, putch
import json
import datetime
import random
import numpy as np

class FluxGenerator:
    def __init__(self, model, startReaction, include, initialExclude):
        self.model = model
        self.reactions = model.reactions
        self.startReaction = startReaction
        self.include = include
        self.exclude = set(initialExclude)
        self._useManual = False
        self._verboseOutput = False
        self._removeDuplicates = False
        self._maxCount = -1
        self._maxTime = -1
        self._extra = lambda flux, excludeSet, **kwargs: None
        self._infeasibleExtra = lambda newest, excludeSet, **kwargs: None
        self._autoStopRatio = -1
        self._autoStopMin = -1
        self._manualStop = True
        self._stopReason = None
        self._countDumpFile = None
        self.alpha = 2.8
        self.beta = 2.9
        self.gamma = 4
        
        # Output
        self.efmsGenerated = [] # List of EFMs generated
        self.infeasibleCount = 0 # Count of infeasible EFMs tried
        self.duplicateCount = 0 # Count of duplicate EFMs generated
        self.timeTaken = 0
        self.solveTime = 0
        
    def getConfig(self):
        """Returns a string showing the configuration of the FluxGenerator"""
        s = """Model: {model!s}
Start reaction: {startReaction!s}
Initial exclude: {exclude!s}
Use manual: {useManual}
Remove duplicates: {removeDuplicates}
Max count: {maxCount}
Max time: {maxTime}
Auto stop ratio: {asr}
Alpha: {alpha}
Beta: {beta}""".format(model = self.model,
                              startReaction = self.startReaction,
                              exclude = self.exclude,
                              useManual = self._useManual,
                              removeDuplicates = self._removeDuplicates,
                              maxCount = self._maxCount,
                              maxTime = self._maxTime,
                              asr = self._autoStopRatio,
                              alpha = self.alpha,
                              beta = self.beta)
        return s
        
        
        
    def useManualInput(self):
        self._useManual = True
        self._verboseOutput = True
        self.output("q to quit")
        self.output("e to view EFMs generated")
        self.output("c to view counts")
        
    def setMaxCount(self, maxCount):
        self.output("Setting max count " + str(maxCount))
        self._maxCount = maxCount
        
    def setMaxTime(self, maxTime):
        self.output("Setting max time " + str(maxTime))
        self._maxTime = maxTime
        
    def useAutoStop(self, ratio=1, minn=8):
        self.output("Using auto stop with ratio=" + str(ratio) + " min=" + str(min))
        self._autoStopRatio = ratio
        self._autoStopMin = minn
        
    def disableManualStop(self):
        self.output("Manual stop turned off")
        self._manualStop = False
        
    def removeDuplicates(self):
        self.output("Removing intermediate duplicates.")
        self._removeDuplicates = True
        
    def setExtra(self, extra):
        self._extra = extra
        
    def setInfeasibleExtra(self, infeasibleExtra):
        self._infeasibleExtra = infeasibleExtra
        
    def verboseOutput(self):
        self._verboseOutput = True
        
    def silent(self):
        self.output = lambda string, **kwargs : None
        
    def stop(self, reason):
        """ Stops the FluxGenerator at the next iteration of the loop."""
        if reason is not None:
            self._stopReason = reason
        self._stopping = True

    
    def getMinimalEFMs(self):
        if len(self.uniqueEFMs) == 0:
            return []
        # Get all bits that correspond to external reactions
        externalReactionIDs = set([r.id for r in self.reactions if r.id[:2] == "EX"])
        
        # Get the unique efms sorted according to their value
        # Remove external reactions
        efms = sorted([fluxToNames(e) - externalReactionIDs for e in self.uniqueEFMs])
        
        minimalEFMs = [efms[0]]
        
                
        # Go through the efms and find minimal EFMs
        for e in efms[1:]:
            minimal = True
            # Remove all bits corresponding to external reactions
            e = e - externalReactionIDs
            for m in minimalEFMs:
                # Try and OR in all previously found minimal EFMs
                orv = e | m
                if orv == e:
                    # This is not minimal
                    minimal = False
                    break
            if minimal:
                minimalEFMs.append(e)
                
        return minimalEFMs
        
    def getTimeDelta(self):
        now = time.time()
        return now - self._startTime - self._waitTime
    
    def _feedback(self, string):
        if self._verboseOutput:
            print(string)
            
    def output(self, string, **kwargs):
        print(string, **kwargs)
            
    def getResults(self):
        s = """
Reason to stop: {reason}

Generated:
    {total:>6} total EFMs
    {unique:>6} unique EFMs
    {dupl:>6} duplicate EFMs
    {infs:>6} infeasible EFMs
    Time taken: {time}""".format(total = len(self.efmsGenerated),
                                 unique = len(self.uniqueEFMs),
                                 dupl = self.duplicateCount,
                                 infs = self.infeasibleCount,
                                 time = self.timeTaken,
                                 reason = self._stopReason)
        return s
            
    def printResults(self):
        self.output(self.getResults())
        
    def printProgress(self):
        if (not self._verboseOutput and not self._useManual):
            s = "Time: {:6.2f} Count: {:7} Infeasible: {:7}".format(time.time() - self._startTime,
                                                                    len(self.efmsGenerated),
                                                                    self.infeasibleCount)
            print(s, end="\r")
            
    def pickNext(self, efmsToExplore, reacScores):
        """Given a list of (nextReaction, exclude) tuples, picks the next reaction
        
        Reassign this to change the search strategy. By default, picks the highest scoring reaction
        
        Return: index of the next tuple to use"""
        def getScore(nextReaction):
            if nextReaction in reacScores:
                return reacScores[nextReaction]
            else:
                return 20
        nextReacs = [efm[0] for efm in efmsToExplore]
        max_i = max(range(len(nextReacs)), key = lambda k: getScore(nextReacs[k]))
        return max_i
    
    def score(self, nextReaction, reacCounts):
        """Calculates the score of a reaction given the reacCounts
        
        reacCounts should be a dictionary which contains unique, feasible, and infeasible counts
        for each reaction."""
        if nextReaction not in reacCounts:
            return 1
        uc = reacCounts[nextReaction]["u"]
        fc = reacCounts[nextReaction]["f"]
        dc = fc - uc # Duplicate count
        ic = reacCounts[nextReaction]["i"]
        
        return float(uc+1)**self.alpha / (dc*self.gamma +ic+1)**self.beta
        
    def dumpCountsToFile(self, fileName=None):
        """Causes the FluxGenerator to dump counts to the file.
        
        Returns the name of the file so it can be opened later"""
        if fileName is None:
            date = datetime.datetime.now().strftime("%Y-%m-%d %H%M%S")
            fileName = date + " {model!s} {startReaction!s} counts.csv".format(model = self.model,
                                                                    startReaction = self.startReaction)
        path = "dumps/" + fileName
        self._countDumpFile = open(path, "w")
        self.output("Dumping to file " + path)
        # Write headers
        self._countDumpFile.write("time,total,infeasible,duplicate,unique\n")
        
        return (path)
            
    def genAll(self, num_slaves=1):
        self.output("Finding EFMs. Push ESC to quit.")
        genThread = threading.Thread(target = self._genAll, args=(num_slaves,))
        genThread.start()
        a = None
        # Break loop if the escape key is pressed
        if self._manualStop:
            while genThread.isAlive():
                a = getch()
                if ord(a) == 27: # Escape key
                    self.stop("Manual stop")
                    break
        
        genThread.join()
        print("                                                      ", end="\r")
        
    def _genAll(self, num_slaves):
        class Slave:
            def __init__(self, model, include):
                self.result = None
                self.thread = None
                self.model = model
                self.include = include
                    
            def run(self, nextReaction, excludeSet):
                self.result = None
                def findEFM(excludeSet):
                    flux = models.findEFM(self.model, self.include, excludeSet, 0)
                    # Return as a tuple to distinguish between not having a result and
                    # not finding an answer
                    self.result = (nextReaction, excludeSet, flux)
                    
                self.thread = threading.Thread(target = findEFM, args=(excludeSet,))
                self.thread.start()
                
            def isAlive(self):
                return self.thread is not None and self.thread.isAlive()
            
            def hasResult(self):
                return self.result is not None
        
        def fluxToNames(flux):
            return set(flux.keys())
        
        def checkAutoStop():
            vc = len(self.uniqueEFMs)
            if self._autoStopRatio > 0 and vc > self._autoStopMin:
                ic = self.infeasibleCount                
                if ic > self._autoStopRatio * vc:
                    self.stop("Auto stop: infeasible-to-unique ratio exceeds limit")
                
        def runExtra(flux, excludeVal, **kwargs):
            startWait = time.time()
            self._extra(flux, excludeVal, **kwargs)
            self._waitTime += time.time() - startWait
            
        def runInfeasibleExtra(newest, excludeSet, **kwargs):
            startWait = time.time()
            self._infeasibleExtra(newest, excludeSet, **kwargs)
            self._waitTime += time.time() - startWait
    
        if self._countDumpFile is not None and self._countDumpFile.closed:
            self.output("Count dump file unavailable")
            self._countDumpFile = None
            
        # Initialise slaves
        slaves = []
        for _ in range(num_slaves):
            slaves.append(Slave(self.model, self.include))
            
        # Set up initial values
        s0 = slaves[0]
        s0.run(self.include.keys()[0], self.exclude)
        # Wait for slave to finish
        while s0.isAlive():
            pass
        flux = s0.result[2]
        self.efmsGenerated = [] # List of EFMs generated
        self.uniqueEFMs = [] # List of unique EFMs
        self.infeasibleCount = 0 # Count of infeasible EFMs tried
        self.duplicateCount = 0 # Count of duplicate EFMs generated
        self.timeTaken = 0
        self._startTime = time.time()
        self._waitTime = 0
        self._stopping = False
        
        # efmsToExplore consists of
        # (reactionToAdd, score, prevExcludeSet)
        # tuples
        efmsToExplore = [] # List of EFMs left to check
        efmsGeneratedNames = [] # List of EFMs generated by reaction names
        reacCounts = {}
        reacScores = {}
        
        def thisScore(r):
            return self.score(r, reacCounts)
        
        # Initialise
        if flux is not None:
            self.efmsGenerated.append(flux)
            fluxNames = fluxToNames(flux)
            efmsGeneratedNames.append(fluxNames)
            # Do anything extra that has been specified by the user
            runExtra(flux , self.exclude, newest=self.startReaction)
            knockOutAble = [r for r in fluxNames if r not in self.include and r not in self.exclude]
            efmsToExplore = [(r, self.exclude) for r in knockOutAble]
            
        activeSlaves = list(slaves)
        asi = 0 # active slave index
        testedExclusions = [self.exclude]
        while len(activeSlaves) > 0:
            # Check criteria to finish early
            if len(efmsToExplore) == 0: # No EFMs left to explore
                self.stop("Explored all EFMs")
            if self._maxCount > 0 and len(self.efmsGenerated) >= self._maxCount: # Reached max count
                self.stop("Generated max number of EFMs")
            # Check whether max time has been exceeded
            if self._maxTime > 0 and time.time() - self._startTime >= self._maxTime:
                self.stop("Max time reached")
            checkAutoStop()
            
            nextSlave = activeSlaves[asi]
            result = None
            
            if nextSlave.hasResult():
                # Get the result
                result = nextSlave.result
            if not nextSlave.isAlive():
                if self._stopping:
                    # If the FluxGenerator is stopping then stop this thread
                    activeSlaves.remove(nextSlave)
                else:
                    # Otherwise assign another reaction
                    r = self.pickNext(efmsToExplore, reacScores)
                    nextReac = efmsToExplore[r]
                    efmsToExplore = efmsToExplore[:r] + efmsToExplore[r+1:] # Remove rth element from list
                    nextReaction = nextReac[0]
                    prevExclude = nextReac[1]
                    nextExclude = prevExclude | set([nextReaction])
                    
                    if nextExclude in testedExclusions:
                        continue
                    testedExclusions.append(nextExclude)
                    self._feedback("Processing exclude set " + str(nextExclude))
                    nextSlave.run(nextReaction, nextExclude)
            
            # Check for user input to pause or give feedback
            if self._useManual:
                startWait = time.time()
                userInput = raw_input()
                if userInput == "q":
                    self.stop("User input")
                    break
                if userInput == "e":
                    self.output(self.efmsGenerated)
                if userInput == "c":
                    self.output("Generated" + str(len(self.efmsGenerated)))
                    self.output("Generated" + str(self.duplicateCount) + "duplicate EFMs") 
                    self.output("Generated" + str(self.infeasibleCount) + "infeasible EFMs")
                self._waitTime += time.time() - startWait
                
            # If a result was received
            if result is not None:
                (nextReaction, nextExclude, flux) = result
                if flux is None:
                    # Infeasible model
                    self.infeasibleCount += 1
                    self._feedback("\tInfeasible")
                    runInfeasibleExtra(nextReaction, nextExclude)
                    if nextReaction not in reacCounts:
                        reacCounts[nextReaction] = {}
                        reacCounts[nextReaction]["f"] = 0
                        reacCounts[nextReaction]["u"] = 0
                        reacCounts[nextReaction]["i"] = 1
                    else:
                        reacCounts[nextReaction]["i"] += 1
                    reacScores[nextReaction] = thisScore(nextReaction)
                    # Check auto stop
                    checkAutoStop()
                else:
                    efmReactionsSet = fluxToNames(flux)
                    self._feedback("\tFlux generated " + str(efmReactionsSet))
                    self.efmsGenerated.append(flux)
                    # Feasible model
                    # Check if this has been generated before
                    unique = True
                    if efmReactionsSet in efmsGeneratedNames:
                        self._feedback("\tDuplicate flux")
                        self.duplicateCount += 1
                        unique = False
                    else:
                        # Add to unique list
                        self.uniqueEFMs.append(flux)
                        efmsGeneratedNames.append(efmReactionsSet)
                        
                    # get score
                    if nextReaction not in reacCounts:
                        reacCounts[nextReaction] = {}
                        reacCounts[nextReaction]["f"] = 1
                        reacCounts[nextReaction]["u"] = 0
                        reacCounts[nextReaction]["i"] = 0
                    else:
                        reacCounts[nextReaction]["f"] += 1
                    if unique:
                        reacCounts[nextReaction]["u"] += 1
                    reacScores[nextReaction] = thisScore(nextReaction)
                        
                    # If this generated a new flux
                    if not self._removeDuplicates or unique:
                        # Find all possible reactions that could be knocked out
                        knockOutAble = [r for r in fluxNames if r not in self.include and r not in nextExclude]
                        efmsToExplore += [(r, nextExclude) for r in knockOutAble]
                        
                    # Do anything extra that has been specified by the user
                    runExtra(flux, nextExclude, newest=nextReaction)
                
            if len(activeSlaves) > 0:
                asi = (asi + 1) % len(activeSlaves)
            
            # Dump counts if necessary
            if self._countDumpFile is not None:
                self._countDumpFile.write("{t},{tc},{ic},{dc},{uc}\n".format(t=self.getTimeDelta(),
                                                     tc=len(self.efmsGenerated),
                                                     ic=self.infeasibleCount,
                                                     dc=self.duplicateCount,
                                                     uc=len(self.uniqueEFMs)))
            self.printProgress()
                
        self.timeTaken = self.getTimeDelta()
        self.printResults()
        putch("\n")
        if self._countDumpFile is not None:
            self._countDumpFile.close()
            
        for s in slaves:
            s.thread.join()

        
    def writeResults(self, fileName=None):
        if fileName is None:
            date = datetime.datetime.now().strftime("%Y-%m-%d %H%M%S")
            fileName = date + " {model!s} {startReaction!s} EFMs.json".format(model = self.model,
                                                                    startReaction = self.startReaction)
        data = json.dumps(self.efmsGenerated)
        with open(fileName, "w") as f:
            f.write(data)
    
            