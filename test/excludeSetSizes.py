class Node:
    def __init__ (self, val, nextNode):
        assert val is not None
        self._val = val
        self._nextNode = nextNode
        
    def __str__(self):
        return self._val.__str__()
        
    def getVal(self):
        return self._val
    
    def getNext(self):
        return self._nextNode

class LinkedList:          
    def __init__ (self, head):
        assert head is not None
        self.head = head
        self.len = 1
        
        node = self.head.getNext()
        while node is not None:
            self.len += 1
            node = node.getNext()
        
    def __len__(self):
        return self.len
    
    def __contains__(self, val):
        node = self.head
        while node is not None:
            if node.getVal() == val:
                return True
            node = node.getNext()
        return False
    
    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        s = "["
        node = self.head
        while node is not None:
            s += str(node)
            node = node.getNext()
            if node is not None:
                s += ", "
        s += "]"
        return s
    
    @staticmethod
    def mkNew(val):
        return LinkedList(Node(val, None))
        
    def push(self, val):
        """Adds an item at the head of the list"""
        newHead = Node(val, self.head)
        self.head = newHead
        self.len += 1
    
    def pop(self):
        """Removes first item in the list"""
        assert self.len > 1
        node = self.head
        self.head = self.head.getNext()
        self.len -= 1
        return node.getVal()
    
import cobra
import random
import time
import sys
    
excludeSetsToGen = 1000000
branchFactor = 100

def setWithReactionNames():
    print "Set structure using reaction names"
    if "model" not in locals():
        model = cobra.io.read_sbml_model("MODEL1108160000.xml")
    reactionIDs = [r.id for r in model.reactions]

    excludeSets = []
    genned = []
    
    def genNextExclude(exclude, branchFactor):
        """Take an exclude set and add a random new reaction"""
        rs = random.sample(reactionIDs, branchFactor)
        for r in rs:
            if r not in exclude:
                nexclude = exclude | set([r])
                excludeSets.append(nexclude)
        genned.append(exclude)
            
    initExclude = set(["Ec_biomass_iJO1366_WT_53p95M"])
    excludeSets.append(initExclude)
    
    startTime = time.time()
    while len(excludeSets) < excludeSetsToGen:
        es = random.sample(excludeSets, 1)[0]
        if es in genned:
            # Skip
            continue
        genNextExclude(es, branchFactor)
        
        print "{:10} sets\r".format(len(excludeSets)),
    
    deltaTime = time.time() - startTime
    
    longest = max(excludeSets, key = lambda k: len(k))
    avg = float(sum([len(es) for es in excludeSets])) / len(excludeSets)
    print "Generated {} exclude sets in {:.3f} seconds.".format(len(excludeSets), deltaTime)
    print "Largest: {} reactions, average: {:.2f} reactions".format(len(longest), avg)
    print "Size of exclude set: {} kB".format(sys.getsizeof(excludeSets)/1000)
    
def setWithCobraReactions():
    print "Set structure using cobra reactions"
    if "model" not in locals():
        model = cobra.io.read_sbml_model("MODEL1108160000.xml")
    reactions = model.reactions

    excludeSets = []
    genned = []
    
    def genNextExclude(exclude, branchFactor):
        """Take an exclude set and add a random new reaction"""
        rs = random.sample(reactions, branchFactor)
        for r in rs:
            if r not in exclude:
                nexclude = exclude | set([r])
                excludeSets.append(nexclude)
        genned.append(exclude)
            
    initExclude = set([reactions.get_by_id("Ec_biomass_iJO1366_WT_53p95M")])
    excludeSets.append(initExclude)
    
    startTime = time.time()
    while len(excludeSets) < excludeSetsToGen:
        es = random.sample(excludeSets, 1)[0]
        if es in genned:
            # Skip
            continue
        genNextExclude(es, branchFactor)
        
        print "{:10} sets\r".format(len(excludeSets)),
    
    deltaTime = time.time() - startTime
    
    longest = max(excludeSets, key = lambda k: len(k))
    avg = float(sum([len(es) for es in excludeSets])) / len(excludeSets)
    print "Generated {} exclude sets in {:.3f} seconds.".format(len(excludeSets), deltaTime)
    print "Largest: {} reactions, average: {:.2f} reactions".format(len(longest), avg)
    print "Size of exclude set: {} kB".format(sys.getsizeof(excludeSets)/1000)
    
def linkedListWithNames():
    if "model" not in locals():
        model = cobra.io.read_sbml_model("MODEL1108160000.xml")
    reactionIDs = [r.id for r in model.reactions]

    excludeSets = []
    genned = []
    
    def genNextExclude(exclude, branchFactor):
        """Take an exclude set and add a random new reaction"""
        rs = random.sample(reactionIDs, branchFactor)
        for r in rs:
            if r not in exclude:
                nexclude = LinkedList(exclude.head)
                nexclude.push(r)
                excludeSets.append(nexclude)
        genned.append(exclude)
            
    initExclude = LinkedList.mkNew("Ec_biomass_iJO1366_WT_53p95M")
    excludeSets.append(initExclude)
    
    startTime = time.time()
    while len(excludeSets) < excludeSetsToGen:
        es = random.sample(excludeSets, 1)[0]
        if es in genned:
            # Skip
            continue
        genNextExclude(es, branchFactor)
        
        print "{:10} sets\r".format(len(excludeSets)),
    
    deltaTime = time.time() - startTime
    
    longest = max(excludeSets, key = lambda k: len(k))
    avg = float(sum([len(es) for es in excludeSets])) / len(excludeSets)
    print "Generated {} exclude sets in {:.3f} seconds.".format(len(excludeSets), deltaTime)
    print "Largest: {} reactions, average: {:.2f} reactions".format(len(longest), avg)
    print "Size of exclude set: {} kB".format(sys.getsizeof(excludeSets)/1000)












