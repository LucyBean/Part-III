'''
Created on Feb 16, 2017

@author: Lucy
'''

import threading
import time
import random

class Slave:
    def __init__(self, name):
        self.result = None
        self.thread = None
        self.name = name
        
    def process(self, val):
        self.result = None
        def runTarget(args):
            print "{name} starting work on {args}".format(name=self.name, args=args)
            time.sleep(random.uniform(2, 4)) # Pretend to do work
            self.result = random.randint(1,10)
            
        self.thread = threading.Thread(target = lambda: runTarget(val))
        self.thread.start()      
        
    def isAlive(self):
        return self.thread is not None and self.thread.isAlive()
    
    def hasResult(self):
        return self.result is not None
    
# Create slaves
slaves = []
for i in range(5):
    slaves.append(Slave(i))
    

activeSlaves = list(slaves)
thingsToDo = range(10)
si = 0
# While there are active slaves
while len(activeSlaves) > 0:
    # Check the next slave
    slave = activeSlaves[si]
    if slave.hasResult():
        print "{name} finished with result {result}".format(name=slave.name, result=slave.result)
    if not slave.isAlive():
        if len(thingsToDo) > 0:
            # Give it a thing to do
            thing = thingsToDo[0]
            thingsToDo = thingsToDo[1:]
            slave.process(thing)
        else:
            activeSlaves.remove(slave)
            
    if len(activeSlaves) > 0:
        si = (si + 1) % len(activeSlaves)
    
# Wait for slaves to finish
    
print "All work complete"
