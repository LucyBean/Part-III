'''
Created on Nov 12, 2016

@author: Lucy
'''
import models
import cobra.test

if __name__ == '__main__':
    model = cobra.test.create_test_model("textbook")
    include = ["SUCDi"]
    exclude = ["FRD7"]
    ignore = ["adp_c", "atp_c", "coa_c", "h2o_c", "h_c", "h_e",
                          "nad_c", "nadh_c", "nadp_c", "nadph_c", "pep_c"]
    
    fluxes = models.process(model, include, exclude, ignore)
    
    models.display("e_coli_core.Core metabolism", fluxes)