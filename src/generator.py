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
        self._verboseOutput = True
        self._removeDuplicates = False
        self._maxCount = -1
        self._maxTime = -1
        self._extra = lambda flux, excludeSet, **kwargs: None
        self._infeasibleExtra = lambda newest, excludeSet, **kwargs: None
        self._autoStopRatio = -1
        self._autoStopMin = -1
        self._manualStop = True
        self._stopReason = None
        
        # Output
        self.efmsGenerated = [] # List of EFMs generated
        self.infeasibleCount = 0 # Count of infeasible EFMs tried
        self.duplicateCount = 0 # Count of duplicate EFMs generated
        self.timeTaken = 0
        
    def getConfig(self):
        """Returns a string showing the configuration of the FluxGenerator"""
        s = """Model: {model!s}
Start reaction: {startReaction!s}
Initial exclude: {exclude!s}
Use manual: {useManual}
Remove duplicates: {removeDuplicates}
Max count: {maxCount}
Max time: {maxTime}
Auto stop ratio: {asr}""".format(model = self.model,
                              startReaction = self.startReaction,
                              exclude = self.exclude,
                              useManual = self._useManual,
                              removeDuplicates = self._removeDuplicates,
                              maxCount = self._maxCount,
                              maxTime = self._maxTime,
                              asr = self._autoStopRatio)
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
        
    def suppressOutput(self):
        self._verboseOutput = False
        
    def silent(self):
        self.output = lambda string, **kwargs : None
        
    def stop(self, reason):
        """ Stops the FluxGenerator at the next iteration of the loop."""
        if reason is not None:
            self._stopReason = reason
        self._breakLoop = True

    
    def getMinimalEFMs(self):
        if len(self.uniqueEFMs) == 0:
            return []
        # Get all bits that correspond to external reactions
        externalReactionIDs = set([r.id for r in self.reactions if r.id[:2] == "EX"])
        
        # Get the unique efms sorted according to their value
        # Remove external reactions
        efms = sorted([self._fluxToNames(e) - externalReactionIDs for e in self.uniqueEFMs])
        
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
    
    def _checkAutoStop(self):
        vc = len(self.uniqueEFMs)
        if self._autoStopRatio > 0 and vc > self._autoStopMin:
            ic = self.infeasibleCount
            
            if ic > self._autoStopRatio * vc:
                self.stop("Auto stop: infeasible-to-unique ratio exceeds limit")
        
    def _runExtra(self, flux, excludeVal, **kwargs):
        startWait = time.time()
        self._extra(flux, excludeVal, **kwargs)
        self._waitTime += time.time() - startWait
        
    def _runInfeasibleExtra(self, newest, excludeSet, **kwargs):
        startWait = time.time()
        self._infeasibleExtra(newest, excludeSet, **kwargs)
        self._waitTime += time.time() - startWait
        
    def getTimeDelta(self):
        now = time.time()
        return now - self._startTime - self._waitTime
    
    def _findEFM(self, exclude):
        return models.findEFM(self.model, self.include, exclude, 0)
    
    def _fluxToNames(self, flux):
        return set(flux.keys())
    
    def _feedback(self, string):
        if self._verboseOutput:
            print(string)
            
    def output(self, string, **kwargs):
        print(string, **kwargs)
            
    def getResults(self):
        s = """Reason to stop: {reason}

Generated:
    {total:>6} total EFMs
    {unique:>6} unique EFMs
    {dupl:>6} duplicate EFMs
    {infs:>6} infeasible EFMs
    {min:>6} minimal EFMs
    Time taken: {time}""".format(total = len(self.efmsGenerated),
                                 unique = len(self.uniqueEFMs),
                                 dupl = self.duplicateCount,
                                 infs = self.infeasibleCount,
                                 min = len(self.getMinimalEFMs()),
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
            
    def genAll(self, strategy):
        self.output("Finding EFMs. Push ESC to quit.")
        if strategy == 1:
            self.output("Using strategy 1")
            genThread = threading.Thread(target=self._genAll)
        elif strategy == 2:
            self.output("Using strategy 2")
            genThread = threading.Thread(target=self._genAllRand)
        else:
            return
        genThread.start()
        a = None
        # Break loop if the escape key is pressed
        if self._manualStop:
            while genThread.isAlive():
                a = getch()
                if ord(a) == 27: # Escape key
                    self.stop("Manual stop")
                    break
        else:
            genThread.join()
            
        self.output("\n")
        
    def _genAll(self):
        # Set up initial values
        flux  = self._findEFM(self.exclude)
        exclude = self.exclude
        self.efmsGenerated = [] # List of EFMs generated
        self.uniqueEFMs = [] # List of unique EFMs
        self.infeasibleCount = 0 # Count of infeasible EFMs tried
        self.duplicateCount = 0 # Count of duplicate EFMs generated
        self.timeTaken = 0
        self._startTime = time.time()
        self._waitTime = 0
        self._breakLoop = False
        
        efmsToExplore = [] # List of EFMs left to check
        efmsGeneratedNames = [] # List of EFMs generated by reaction names
        
        if flux  is not None:
            efmsToExplore.append((flux , self.exclude))
            self.efmsGenerated.append(flux )
            efmsGeneratedNames.append(self._fluxToNames(flux ))
            # Do anything extra that has been specified by the user
            self._runExtra(flux , self.exclude, newest=self.startReaction)
            
        testedExclusions = [self.exclude]
        
        while not self._breakLoop:
            if efmsToExplore == []: # No EFMs left to explore
                self.stop("Explored all EFMs")
                break
            if self._maxCount > 0 and len(self.efmsGenerated) >= self._maxCount: # Reached max count
                self.stop("Generated max number of EFMs")
                break
            # Check whether max time has been exceeded
            if self._maxTime > 0 and time.time() - self._startTime >= self._maxTime:
                self.stop("Max time reached")
                break
            
            efm = efmsToExplore[-1] # Get the last element in the list
            efmsToExplore = efmsToExplore[:-1]
            flux  = efm[0] # Extract the flux (list of reactions)
            fluxNames = self._fluxToNames(flux ) # Convert list of reactions to val
            exclude = efm[1] # Get the exclusion list that generated the flux
            self._feedback("Exploring EFM " + str(fluxNames))
            
            # Find all possible knock outs
            knockOutAble = [self.reactions.get_by_id(r) for r in fluxNames
                            if r not in self.include and r not in exclude]
            # And try to knock out each reaction
            for r in knockOutAble:
                nexcludeSet = exclude | set([r.id])
                self._feedback("Next exclusion set " + str(nexcludeSet))
                
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
                    
                # Check if this exclusion set has been tested before
                if nexcludeSet in testedExclusions:
                    self._feedback("Duplicate exclusion set " + str(nexcludeSet))
                    continue # Skip this exclusion set and try next
                testedExclusions.append(nexcludeSet)
                
                # Generate the EFM
                nflux = self._findEFM(nexcludeSet)
                
                if nflux is None:
                    # Infeasible model
                    self.infeasibleCount += 1
                    self._feedback("\tInfeasible")
                    self._infeasibleExtra(r, nexcludeSet)
                    self.printProgress()
                    # Check auto stop
                    self._checkAutoStop()
                    # Check if stop has been flagged
                    if self._breakLoop:
                        break
                    # Check whether max time has been exceeded
                    if self._maxTime > 0 and self.getTimeDelta() >= self._maxTime:
                        break
                else:
                    # Feasible model
                    efmReactionsSet = self._fluxToNames(nflux)
                    self._feedback("\tFlux generated " + str(efmReactionsSet))
                    
                    # Check if this has been generated before
                    unique = True
                    if efmReactionsSet in efmsGeneratedNames:
                        self._feedback("\tDuplicate flux")
                        self.duplicateCount += 1
                        unique = False
                    else:
                        # Add to unique list
                        self.uniqueEFMs.append(nflux)
                        efmsGeneratedNames.append(efmReactionsSet)
                    
                    # Add this to the list of EFMs to check if:
                        # not removing duplicates, or
                        # it is not a duplicate
                    if not self._removeDuplicates or unique:
                        efmsToExplore.append((nflux, nexcludeSet))
                    
                    # Add to list of EFMs generated
                    self.efmsGenerated.append(nflux)
                        
                    # Do anything extra that has been specified by the user
                    self._runExtra(nflux, nexcludeSet, newest=r)
                    self.printProgress()
                    
                    # Check auto stop
                    self._checkAutoStop()
                    # Check if stop has been flagged
                    if self._breakLoop:
                        break
                    
                    # Check whether max number has been generated
                    if self._maxCount > 0 and len(self.efmsGenerated) >= self._maxCount:
                        break
                    
                    # Check whether max time has been exceeded
                    if self._maxTime > 0 and time.time() - self._startTime >= self._maxTime:
                        break
                
        self.timeTaken = self.getTimeDelta()
        self.printResults()
        putch("\n")
        
    def score(self, r, reacCounts):
        if r not in reacCounts:
            return 1
        fc = reacCounts[r]["f"]
        ic = reacCounts[r]["i"]
        return (fc+1)/float(fc+ic+1)
        
    def _genAllRand(self):
        # Set up initial values
        flux  = self._findEFM(self.exclude)
        exclude = self.exclude
        self.efmsGenerated = [] # List of EFMs generated
        self.uniqueEFMs = [] # List of unique EFMs
        self.infeasibleCount = 0 # Count of infeasible EFMs tried
        self.duplicateCount = 0 # Count of duplicate EFMs generated
        self.timeTaken = 0
        self._startTime = time.time()
        self._waitTime = 0
        self._breakLoop = False
        
        # efmsToExplore consists of
        # (reactionToAdd, score, prevExcludeSet)
        # tuples
        
        efmsToExplore = [] # List of EFMs left to check
        efmsGeneratedNames = [] # List of EFMs generated by reaction names
        reacCounts = {}
        
        def thisScore(r):
            return self.score(r, reacCounts)
        
        # Initialise
        if flux  is not None:
            efmsToExplore.append((flux , 1, set(exclude)))
            self.efmsGenerated.append(flux )
            fluxNames = self._fluxToNames(flux)
            efmsGeneratedNames.append(fluxNames)
            # Do anything extra that has been specified by the user
            self._runExtra(flux , exclude, newest=self.startReaction)
            knockOutAble = [r for r in fluxNames if r not in self.include and r not in exclude]
            efmsToExplore = [(r, 1, exclude) for r in knockOutAble]
            
        testedExclusions = [self.exclude]
        
        while not self._breakLoop:
            if len(efmsToExplore) == 0: # No EFMs left to explore
                self.stop("Explored all EFMs")
                break
            if self._maxCount > 0 and len(self.efmsGenerated) >= self._maxCount: # Reached max count
                self.stop("Generated max number of EFMs")
                break
            # Check whether max time has been exceeded
            if self._maxTime > 0 and time.time() - self._startTime >= self._maxTime:
                self.stop("Max time reached")
                break
            self._checkAutoStop()
            if self._breakLoop:
                break
            
            scores = np.array([e[1] for e in efmsToExplore])
            totScore = float(np.sum(scores))
            relScores = scores / totScore
            rand = random.uniform(0,1)
            r = np.sum(np.cumsum(relScores) < rand)
            nextReac = efmsToExplore[r]
            efmsToExplore = efmsToExplore[:r] + efmsToExplore[r+1:] # Remove rth element from list
            nextReaction = nextReac[0]
            prevExclude = nextReac[2]
            nextExclude = prevExclude | set([nextReaction])
            
            if nextExclude in testedExclusions:
                self._feedback("Duplicate exclusion set " + str(nextExclude))
                continue # Skip this exclusion set and try next
            testedExclusions.append(nextExclude)
                
            self._feedback("Trying exclude set " + str(nextExclude))
            flux = self._findEFM(nextExclude)
            self._feedback("Success: " + str(flux is not None))
            
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
                
            if flux is None:
                # Infeasible model
                self.infeasibleCount += 1
                self._feedback("\tInfeasible")
                self._infeasibleExtra(nextReac, nextExclude)
                if nextReaction  not in reacCounts:
                    reacCounts[nextReaction] = {}
                    reacCounts[nextReaction]["f"] = 0
                    reacCounts[nextReaction]["u"] = 0
                    reacCounts[nextReaction]["i"] = 1
                else:
                    reacCounts[nextReaction]["i"] += 1
                # Check auto stop
                self._checkAutoStop()
            else:
                efmReactionsSet = self._fluxToNames(flux)
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
                    
                # If this generated a new flux
                if not self._removeDuplicates or unique:
                    # Find all possible reactions that could be knocked out
                    knockOutAble = [r for r in fluxNames if r not in self.include and r not in exclude]
                    efmsToExplore += [(r, thisScore(r), nextExclude) for r in knockOutAble]
                    
                # Do anything extra that has been specified by the user
                self._runExtra(flux, nextExclude, newest=r)
                
            self.printProgress()
                
        self.timeTaken = self.getTimeDelta()
        self.printResults()
        putch("\n")
        
    def writeResults(self, label=""):
        date = datetime.datetime.now().strftime("%Y-%m-%d %H%M%S")
        reac = self.startReaction.id
        data = json.dumps(self.efmsGenerated)
        with open(date + " " + label + " " + reac + ".json", "w") as f:
            f.write(data)
    
            