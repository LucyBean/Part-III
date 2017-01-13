'''
Created on Jan 12, 2017

@author: Lucy
'''

from src import models, helpers
import time

class FluxGenerator:
    def __init__(self, model, startReaction, include, initialExclude):
        self.model = model
        self.reactions = model.reactions
        self.reactionVals = helpers.genReactionVals(self.reactions)
        self.startReaction = startReaction
        self.include = include
        self.exclude = initialExclude
        self._useManual = False
        self._verboseOutput = True
        self._removeDuplicates = False
        self._maxCount = -1
        self._extra = lambda flux, excludeVal: None
        
        # Output
        self.efmsGenerated = [] # List of EFMs generated
        self.infeasibleCount = 0 # Count of infeasible EFMs tried
        self.duplicateCount = 0 # Count of duplicate EFMs generated
        self.timeTaken = 0
        
    def useManualInput(self):
        self._useManual = True
        self._verboseOutput = True
        print "q to quit"
        print "e to view EFMs generated"
        print "c to view counts"
        
    def setMaxCount(self, maxCount):
        self._maxCount = maxCount
        
    def removeDuplicates(self):
        self._removeDuplicates = True
        
    def setExtra(self, extra):
        self._extra = extra
        
    def suppressOutput(self):
        self._verboseOutput = False
        
    def _reactionNamesToVal(self, reactionIDs):
        return helpers.reactionNamesToVal(reactionIDs, self.reactionVals)
    
    def _findEFM(self, exclude):
        return models.findEFM(self.model, self.include, exclude, 0)
    
    def _output(self, string):
        if self._verboseOutput:
            print string
            
    def printResults(self):
        print "Generated:"
        print "\t", len(self.efmsGenerated), "total EFMs"
        print "\t", self.duplicateCount, "duplicate EFMs generated"
        print "\t", self.infeasibleCount, "infeasible EFMs tried"
        print "Time taken:", self.timeTaken
                      
    def genAll(self):
        # Set up initial values
        flux = self._findEFM(self.exclude)
        self.efmsGenerated = [] # List of EFMs generated
        self.infeasibleCount = 0 # Count of infeasible EFMs tried
        self.duplicateCount = 0 # Count of duplicate EFMs generated
        self.timeTaken = 0
        
        efmsToExplore = [] # List of EFMs left to check
        
        if flux is not None:
            efmsToExplore.append((flux, self.exclude))
            self.efmsGenerated.append(self._reactionNamesToVal(flux))
            # Do anything extra that has been specified by the user
            self._extra(flux, self._reactionNamesToVal(self.exclude))
        
        startTime = time.time()
        waitTime = 0
        userInput = ""
        testedExclusions = [self._reactionNamesToVal(self.exclude)]
        
        while userInput != "q":
            if efmsToExplore == []: # No EFMs left to explore
                self._output("\nExplored all EFMs")
                break
            if self._maxCount > 0 and len(self.efmsGenerated) >= self._maxCount: # Reached max count
                self._output("\nGenerated max number of EFMs")
                break
            
            efm = efmsToExplore[-1] # Get the last element in the list
            efmsToExplore = efmsToExplore[:-1]
            flux = efm[0] # Extract the flux (list of reactions)
            fluxVal = self._reactionNamesToVal(flux) # Convert list of reactions to val
            exclude = efm[1] # Get the exclusion list that generated the flux
            self._output("Exploring EFM " + str(fluxVal))
            
            # Find all possible knock outs
            knockOutAble = [r for r in self.reactions if r.id not in self.include and r.id not in exclude]
            # And try to knock out each reaction
            for r in knockOutAble:
                nexclude = exclude + [r.id]
                nexcludeVal = self._reactionNamesToVal(nexclude)
                self._output("Next exclusion set " + str(nexcludeVal))
                
                # Check for user input to pause or give feedback
                if self._useManual:
                    startWait = time.time()
                    userInput = raw_input()
                    if userInput == "q":
                        break
                    if userInput == "e":
                        print self.efmsGenerated
                    if userInput == "c":
                        print "Generated", len(self.efmsGenerated)
                        print "Generated", self.duplicateCount, "duplicate EFMs" 
                        print "Generated", self.infeasibleCount, "infeasible EFMs"
                    waitTime += time.time() - startWait
                    
                # Check if this exclusion set has been tested before
                if nexcludeVal in testedExclusions:
                    self._output("Duplicate exclusion set " + str(nexcludeVal))
                    continue # Skip this exclusion set and try next
                testedExclusions.append(nexcludeVal)
                
                # Generate the EFM
                nflux = self._findEFM(nexclude)
                
                if nflux is None:
                    # Infeasible model
                    self.infeasibleCount += 1
                    self._output("\tInfeasible")
                else:
                    # Feasible model
                    efmReactions = nflux.keys()
                    efmReactionsVal = self._reactionNamesToVal(efmReactions)
                    self._output("\tFlux generated " + str(efmReactionsVal))
                    
                    # Check if this has been generated before
                    if efmReactionsVal in self.efmsGenerated:
                        self.duplicateCount += 1
                    
                    # Add this to the list of EFMs to check if:
                        # not removing duplicates, or
                        # it is not a duplicate
                    if not self._removeDuplicates or efmReactionsVal not in self.efmsGenerated:
                        self._output("\tFound flux " + str(efmReactionsVal))
                        efmsToExplore.append((nflux, nexclude))
                        self.efmsGenerated.append(efmReactionsVal)
                        
                    # Do anything extra that has been specified by the user
                    self._extra(nflux, nexcludeVal)
                    
                    # Check whether max number has been generated
                    if self._maxCount > 0 and len(self.efmsGenerated) >= self._maxCount:
                        break
                
        endTime = time.time()
        self.timeTaken = endTime - startTime - waitTime
        self.printResults()
        
    def getMinimalEFMs(self):
        # Get the unique efms sorted according to their value
        efms = sorted(set(self.efmsGenerated))
        minimalEFMs = [efms[0]]
        
        # Go through the efms and find minimal EFMs
        for e in efms[1:]:
            minimal = True
            for m in minimalEFMs:
                # Try and OR in all previously found minimal EFMs
                orv = e | m
                if orv == e:
                    # This is not minimal
                    minimal = False
                    break
            if minimal:
                minimalEFMs.append(e)
                
        return minimal
            
            
            
            
            
            