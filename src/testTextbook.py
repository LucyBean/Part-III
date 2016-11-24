'''
Created on Nov 12, 2016

@author: Lucy
'''
import models
import display
import cobra.test

if __name__ == '__main__':
    model = cobra.test.create_test_model("textbook")
    include = {"EX_lac__D_e":models.FORWARD}
    exclude = []
    #ignore = ["adp_c", "atp_c", "coa_c", "h2o_c", "h_c", "h_e", "nad_c", "nadh_c", "nadp_c", "nadph_c"]
    ignore = []
    
    fluxes = models.findEFM(model, include, exclude, ignore)
    
    if fluxes is not None:
        display.displayEFM(map_name="e_coli_core.Core metabolism", reaction_data=fluxes)