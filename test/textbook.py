'''
Created on Nov 12, 2016

@author: Lucy
'''
from src import models, display
import cobra.test

if __name__ == '__main__':
    model = cobra.test.create_test_model("textbook")
    include = {"EX_glc__D_e":models.REVERSE}
    exclude = ["PGI"]
    
    fluxes = models.findEFM(model, include, exclude)
    
    if fluxes is not None:
        display.displayEFM(map_name="e_coli_core.Core metabolism", reaction_data=fluxes)