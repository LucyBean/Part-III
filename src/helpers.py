'''
Created on Nov 25, 2016

@author: Lucy
'''

def printFlux(flux):
    sortedKeys = sorted(flux)
    for f in sortedKeys:
        print "\t", f, "\t%.3f" % flux[f]
        
def printProducts(products):
    for p in products:
        print p
        for f in products[p]:
            printFlux(f)
            print ""
        print "\n-----\n"