'''
Created on Nov 15, 2016

@author: Lucy
'''

import cobra
from src import models, display

possibleTerminals = ["G6P", "R5Pex", "Pyr"]
startID = "G6P"
    
cobraModel = cobra.io.load_json_model("toyModel.json")

products = models.findProducts(cobraModel, startID, possibleTerminals)
title = "Possible products for starting metabolite " + startID
display.displayAll("toyModelMap.json", products, title)