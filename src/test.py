'''
Created on Nov 12, 2016

@author: Lucy
'''
import models
import cobra.test

if __name__ == '__main__':
    model = cobra.test.create_test_model("textbook")
    include = ["FRUpts2"]
    exclude = []
    ignore = ["adp_c", "atp_c", "coa_c", "h2o_c", "h_c", "h_e",
                          "nad_c", "nadh_c", "nadp_c", "nadph_c", "pep_c"]
    
    fluxes = models.process(model, include, exclude)
    
    if fluxes is not None:
        models.display("e_coli_core.Core metabolism", fluxes)