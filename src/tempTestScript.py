'''
Created on Nov 15, 2016

@author: Lucy
'''

import cobra
import models
import display

externalMetabolites = ["NADP", "NADPH", "CO2", "G6P", "R5Pex", "ATP", "ADP", "Pyr", "NAD", "NADH"]
possibleTerminals = ["G6P", "R5Pex", "Pyr"]
startID = "G6P"
    
cobraModel = cobra.io.load_json_model("toyModel.json")

products = models.findProducts(cobraModel, startID, possibleTerminals)
title = "Possible products for starting metabolite " + startID
display.displayAll(cobraModel, "toyModelMap.json", products, title)