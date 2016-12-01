'''
Created on Nov 25, 2016

@author: Lucy
'''

def printFlux(flux):
    sortedKeys = sorted(flux)
    for f in sortedKeys:
        print "\t", f, "\t%.3f" % flux[f]
        
def printProducts(products):
    for name in products:
        print name
        fluxes = products[name]["fluxes"]
        for flux in fluxes:
            printFlux(flux)
            print ""
        print "\n-----\n"