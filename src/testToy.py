'''
Created on Nov 14, 2016

@author: Lucy
'''

import models
import cobra

defaultIgnore = ["NADP", "NADPH", "CO2", "G6P", "R5Pex", "ATP", "ADP", "Pyr", "NAD", "NADH"]

def testEFM(include, exclude = [], ignore=defaultIgnore):
    model = cobra.io.load_json_model("toyModel.json")
    
    fluxes = models.findEFM(model, include, exclude, ignore)
    if fluxes is not None:
        models.displayEFM(map_json="toyModelMap.json", reaction_data=fluxes)
        
def testPsemi(include):
    model = cobra.io.load_json_model("toyModel.json")
    
    fluxes = models.findPsemi(model, include)
    if fluxes is not None:
        models.displayPsemi(map_json="toyModelMap.json", metabolite_data=fluxes)
        
def emfA():
    include = {"Pgi": models.FORWARD, "Pyk": models.FORWARD}
    exclude = ["Prs_DeoB", "Zwf"]
    testEFM(include, exclude)
    
def emfB():
    include = {"Zwf": models.FORWARD}
    exclude = ["Prs_DeoB", "Pfk", "Fbp"]
    testEFM(include, exclude)
    
def emfC():
    include = {"Zwf": models.FORWARD}
    exclude = ["Prs_DeoB", "Pgi"]
    testEFM(include, exclude)
    
def emfD():
    include = {"Zwf": models.FORWARD}
    exclude = []
    testEFM(include, exclude)
    
def emfE():
    include = {"Pgi": models.FORWARD}
    exclude = []
    testEFM(include, exclude)
    
def emfF():
    include = {"Zwf": models.FORWARD}
    exclude = ["Prs_DeoB"]
    testEFM(include, exclude)
    
def emfG():
    include = {"Pfk": models.FORWARD}
    exclude = []
    testEFM(include, exclude)
    
def ru5p_semi():
    include = ["Ru5P"]
    testPsemi(include)

if __name__ == '__main__':
    ru5p_semi()