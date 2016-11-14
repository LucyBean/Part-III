'''
Created on Nov 14, 2016

@author: Lucy
'''

import models
import cobra

if __name__ == '__main__':
    model = cobra.io.load_json_model("toyModel.json")
    include = {"Zwf":models.FORWARD}
    exclude = ["Prs_DeoB", "Pfk", "Fbp"]
    ignore = ["NADP", "NADPH", "CO2", "G6P", "R5Pex", "ATP", "ADP", "Pyr", "NAD", "NADH"]
    
    fluxes = models.process(model, include, exclude, ignore)
    if fluxes is not None:
        models.display(map_json="toyModelMap.json", reaction_data=fluxes)

