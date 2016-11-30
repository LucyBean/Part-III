'''
Created on Nov 15, 2016

@author: Lucy
'''

import cobra
from src import models, display

startID = "G6P"
    
cobraModel = cobra.io.load_json_model("toyModel.json")

products = models.findProducts(cobraModel, startID)
title = "Possible products for starting metabolite " + startID
display.displayAll(map_json="toyModelMap.json", toDisplay=products, title=title)