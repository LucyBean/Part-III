'''
Created on Nov 15, 2016

@author: Lucy
'''

import cobra
import sys
import models

if __name__ == '__main__':
    cobraModel = cobra.io.load_json_model("toyModel.json")
    
    # EFM A
    reactionsToInclude = {"Pgi": models.FORWARD, "Pyk": models.FORWARD}
    reactionsToExclude = ["Prs_DeoB", "Zwf"]
    
    #EFM D
    #reactionsToInclude = {"Zwf": models.FORWARD}
    #reactionsToExclude = []
    
    #EFM F
    #reactionsToInclude = {"Zwf": models.FORWARD}
    #reactionsToExclude = ["Prs_DeoB"]
    
    externalMetabolites = ["NADP", "NADPH", "CO2", "G6P", "R5Pex", "ATP", "ADP", "Pyr", "NAD", "NADH"]
    
    
    flux = models.findEFM(cobraModel, reactionsToInclude, reactionsToExclude, externalMetabolites)
    reactions = cobraModel.reactions

    possibleTerminals = ["G6P", "R5Pex", "Pyr"]
    (terminalReactants, terminalProducts) = models.findTerminalReactantsAndProducts(possibleTerminals,
                                                                             flux, reactions)

    print "Products:", terminalProducts
    print "Reactants:", terminalReactants
