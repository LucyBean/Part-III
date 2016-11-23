'''
Created on Nov 15, 2016

@author: Lucy
'''

import sys
import cobra
import models


externalMetabolites = ["NADP", "NADPH", "CO2", "G6P", "R5Pex", "ATP", "ADP", "Pyr", "NAD", "NADH"]
possibleTerminals = ["G6P", "R5Pex", "Pyr"]
startID = "G6P"
    
cobraModel = cobra.io.load_json_model("toyModel.json")
metabolites = cobraModel.metabolites
reactions = cobraModel.reactions
if startID not in metabolites:
    sys.stderr.write("Attempting to start from unknown metabolite " + startID)
    return
startingMetabolite = metabolites.get_by_id(startID)
possibleStartReactions = list(startingMetabolite.reactions)
products = {}

for startReaction in possibleStartReactions:
    reactionsToInclude = {startReaction.id: models.FORWARD}
    reactionsToExclude = []
    
    # repeat 10 times
    for i in range(10):
        flux = models.findEFM(cobraModel, reactionsToInclude, reactionsToExclude, externalMetabolites)
        
        if flux is None:
            print "Flux is none."
            break
        
        
        (terminalReactants, terminalProducts) = models.findTerminalReactantsAndProducts(possibleTerminals, flux, reactions)
        
        # Take the first product in terminalProducts and knock out that reaction
        if terminalProducts != []:
            terminalProduct = metabolites.get_by_id(terminalProducts[0])
            rs = list(terminalProduct.reactions)
            # Take the first reaction in this list
            knockOut = [r for r in rs if r.id in flux][0]
            reactionsToExclude.append(knockOut.id)
        # The EMF produced is a loop
        else:
            terminalProduct = startingMetabolite
            rs = list(terminalProduct.reactions)
            # Find all the reactions that the starting metabolite is involved in
            mightKnockOut = [r for r in rs if r.id in flux]
            # Find all of the ones that produce the starting metabolite
            knockOut = []
            for r in mightKnockOut:
                if flux[r.id] > 0 and startingMetabolite in r.products:
                    knockOut.append(r)
                elif flux[r.id] < 0 and startingMetabolite in r.reactants:
                    knockOut.append(r)
            reactionsToExclude.append(knockOut[0].id)
        
        # Add this possible product/pathway to the products dict
        tpid = terminalProduct.id
        if tpid not in products:
            products[tpid] = []
        products[tpid].append(flux)
        
for p in products:
    print "\n", p
    for f in products[p]:
        print "\t", f

