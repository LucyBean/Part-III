'''
Created on Nov 15, 2016

@author: Lucy
'''

# hack hack hack
if "products" not in locals():
    execfile("tempTestScript.py")

possibleTerminals = ["G6P", "R5Pex", "Pyr"]

# extract the metabolites that are only produced or consumed
# use sets to prevent extra work when handling duplicates
onlyProduced = set([m for m in products if m not in reactants])
onlyConsumed = set([m for m in reactants if m not in products])

terminalReactants = [m for m in onlyConsumed if m in possibleTerminals]
terminalProducts =  [m for m in onlyProduced if m in possibleTerminals]
    