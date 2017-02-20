'''
Created on Jan 12, 2017

@author: Lucy
'''

from __future__ import print_function
from src import models
import time
from msvcrt import getch, putch
import json
import datetime
import threading
import cobra

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
        
        self._reset()
        
    def _reset(self):
         # Output
        self.infeasibleCount = 0 # Count of infeasible EFMs tried
        self.duplicateCount = 0 # Count of duplicate EFMs generated
        self.timeTaken = 0
        self.solveTime = 0
        
        # Things needed for solving
        # toTry consists of {nextReaction : [nexcludeSets]}
        self.toTry = {} # Lists of excludeSets to try indexed by next reaction
        self.toTryNum = 0 # Number of things to try
        self.resultsByNames = [] # List of EFMs generated by reaction names
        self.reacCounts = {}
        self.reacScores = {}
        self.testedExclusions = []
        self.uniqueEFMs = [] # List of unique EFMs
        self.infeasibleCount = 0 # Count of infeasible EFMs tried
        self.duplicateCount = 0 # Count of duplicate EFMs generated
        self.timeTaken = 0
        self._startTime = time.time()
        self._waitTime = 0
        self._stopping = False
        
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
            
    def output(self, string, **kwargs):
        print(string, **kwargs)
            
    def getResults(self):
        s = """
Reason to stop: {reason}

Generated:
    {unique:>6} unique EFMs
    {dupl:>6} duplicate EFMs
    {infs:>6} infeasible EFMs
    Time taken: {time}""".format(unique = len(self.uniqueEFMs),
                                 dupl = self.duplicateCount,
                                 infs = self.infeasibleCount,
                                 time = self.timeTaken,
                                 reason = self._stopReason)
        return s
            
    def printResults(self):
        self.output(self.getResults())
        
    def printProgress(self):
        if (not self._verboseOutput and not self._useManual):
            s = ("Time: {:6.2f} Unique: {:7} Duplicate: {:7} Infeasible: {:7}"
                    .format(time.time() - self._startTime,
                    len(self.uniqueEFMs),
                    self.duplicateCount,
                    self.infeasibleCount))
            print(s, end="\r")
            
    def pickNext(self):
        """Returns the next tuple to try.
        Splices it out of the toTry list"""
        def getScore(nextReaction):
            if nextReaction in self.reacScores:
                return self.reacScores[nextReaction]
            else:
                # Arbitrary default
                return 20
        
        # Sort all the next possible reactions by score
        nextReacs = sorted(self.toTry.keys(), key = lambda k: getScore(k))
        # Only keep those that actually have exclude sets left to try
        nextReacs = [nr for nr in nextReacs if len(self.toTry[nr]) > 0]
        i = 0
        # Find the next exclude set for this reaction that has not been tried
        while self.toTryNum > 0:
            while len(self.toTry[nextReacs[i]]) == 0:
                # Find first non empty
                i += 1
                if i == len(self.toTry):
                    break
            if i == len(self.toTry):
                assert self.toTryNum == 0
                break
            excludeSets = self.toTry[nextReacs[i]]
            es = excludeSets[0]
            excludeSets = self.toTry[nextReacs[i]] = excludeSets[1:]
            self.toTryNum -= 1
            if es not in self.testedExclusions:
                return nextReacs[i], es
            
                        
    def addToTryList(self, reaction, exclude):
        if reaction not in self.toTry:
            self.toTry[reaction] = []
        self.toTry[reaction].append(exclude)
        self.toTryNum += 1
    
    def score(self, nextReaction):
        """Calculates the score of a reaction given the reacCounts
        
        reacCounts should be a dictionary which contains unique, feasible, and infeasible counts
        for each reaction."""
        if nextReaction not in self.reacCounts:
            return 1
        uc = self.reacCounts[nextReaction]["u"]
        dc = self.reacCounts[nextReaction]["d"]
        ic = self.reacCounts[nextReaction]["i"]
        
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
            
    def genAll(self):
        self.output("Finding EFMs. Push ESC to quit.")
        genThread = threading.Thread(target = self._genAll)
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
        
    def _genAll(self):
            
        def feedback(string):
            if self._verboseOutput:
                print(string)
        
        def findSolution(excludeSet):
            flux = models.findEFM(self.model, self.include, excludeSet, 0)
            return flux
        
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
        
        def addCount(name, count):
            if name not in self.reacCounts:
                self.reacCounts[name] = {}
                self.reacCounts[name]["d"] = 0
                self.reacCounts[name]["u"] = 0
                self.reacCounts[name]["i"] = 0
            for k in count: 
                self.reacCounts[name][k] += count[k]
                
        
        def genNextExcludeSets(exclude, flux):
            fluxNames = fluxToNames(flux)
            knockOutAble = [r for r in fluxNames if r not in self.include and r not in exclude]
            for r in knockOutAble:
                nexclude = exclude | set([r])
                self.addToTryList(r, nexclude)
                
        def processResult(flux, nextReaction, exclude):
            if flux is None:
                # No solution found
                self.infeasibleCount += 1
                runInfeasibleExtra(nextReaction, exclude)
                feedback("\tInfeasible")
                addCount(nextReaction, {"i":1})
                self.reacScores[nextReaction] = self.score(nextReaction)
                return
            
            names = fluxToNames(flux)
            feedback("\tFlux generated {}".format(names))
            unique = True
            if names in self.resultsByNames:
                # Duplicate result
                feedback("\tDuplicate")
                addCount(nextReaction, {"d":1})
                self.duplicateCount += 1
                unique = False
            else:
                # New result
                self.uniqueEFMs.append(flux)
                addCount(nextReaction, {"u":1})
                self.resultsByNames.append(names)
                
            # Update score
            self.reacScores[nextReaction] = self.score(nextReaction)
                
            # Process this result
            if unique or not self._removeDuplicates:
                genNextExcludeSets(exclude, flux)
                        
        # Set up variables
        if self._countDumpFile is not None and self._countDumpFile.closed:
            self.output("Count dump file unavailable")
            self._countDumpFile = None
            
        self._reset()
         
        # Initialise
        flux = findSolution(self.exclude)
        processResult(flux, self.include.keys()[0], self.exclude)
        self.testedExclusions = [self.exclude]
        
        while True:
            # Check stopping criteria
            if self.toTryNum == 0:
                self.stop("Explored all EFMs")
            numAttempts = len(self.uniqueEFMs) + self.infeasibleCount + self.duplicateCount
            if self._maxCount > 0 and numAttempts >= self._maxCount:
                self.stop("Reached maximum number of attempts")
            if self._maxTime > 0 and time.time() - self._startTime >= self._maxTime:
                self.stop("Max time reached")
            checkAutoStop()
            if self._stopping:
                break
            
            # Pick the next reaction to try
            toTry = self.pickNext()
            assert toTry is not None
            (nextReaction, nexclude) = toTry
            flux = findSolution(nexclude)
            processResult(flux, nextReaction, nexclude)
            
            # Check for manual input
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
        
    def writeResults(self, fileName=None):
        if fileName is None:
            date = datetime.datetime.now().strftime("%Y-%m-%d %H%M%S")
            fileName = date + " {model!s} {startReaction!s} EFMs.json".format(model = self.model,
                                                                    startReaction = self.startReaction)
        data = json.dumps(self.efmsGenerated)
        with open(fileName, "w") as f:
            f.write(data)