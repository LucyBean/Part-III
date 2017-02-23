import cobra
import random

class Stash:
    def __init__(self, reaction, excludeSet):
        self.reaction = reaction
        self.excludeSet = excludeSet
        
    def __eq__(self, x):
        if type(x) == set:
            return self.excludeSet == x
        else:
            return False
        
    def __str__(self):
        s = str(self.reaction) + " " + str(self.excludeSet)
        return s
        
if "model" not in locals():
    model = cobra.io.read_sbml_model("MODEL1108160000.xml")
reactionIDs = [r.id for r in model.reactions]

esets = []
stashes = []
for _ in range(10):
    eset = set(random.sample(reactionIDs, 10))
    sr = random.sample(eset, 1)[0]
    stash = Stash(sr, eset)
    
    stashes.append(stash)
    esets.append(eset)
    
for _ in range(10):
    eset = set(random.sample(reactionIDs, 10))
    sr = random.sample(eset, 1)[0]
    esets.append(eset)  
    
random.shuffle(esets)
    

