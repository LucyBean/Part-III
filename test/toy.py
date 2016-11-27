'''
Created on Nov 14, 2016

@author: Lucy
'''

from src import models, display
import cobra

defaultIgnore = ["NADP", "NADPH", "CO2", "G6P", "R5Pex", "ATP", "ADP", "Pyr", "NAD", "NADH"]
model = cobra.io.load_json_model("toyModel.json")
        
def testPsemi(include):
    model = cobra.io.load_json_model("toyModel.json")
    
    fluxes = models.findPsemi(model, include)
    if fluxes is not None:
        display.displayPsemi(map_json="toyModelMap.json", metabolite_data=fluxes)
        
def emfA():
    include = {"Pgi": models.FORWARD, "Pyk": models.FORWARD}
    exclude = ["Prs_DeoB", "Zwf"]
    return models.findEFM(model, include, exclude)
    
def emfB():
    include = {"Zwf": models.FORWARD}
    exclude = ["Prs_DeoB", "Pfk", "Fbp"]
    return models.findEFM(model, include, exclude)
    
def emfC():
    include = {"Zwf": models.FORWARD}
    exclude = ["Prs_DeoB", "Pgi"]
    return models.findEFM(model, include, exclude)
    
def emfD():
    include = {"Zwf": models.FORWARD}
    exclude = ["Pgi"]
    return models.findEFM(model, include, exclude)
    
def emfE():
    include = {"Pgi": models.FORWARD}
    exclude = []
    return models.findEFM(model, include, exclude)
    
def emfF():
    include = {"Zwf": models.FORWARD}
    exclude = ["Prs_DeoB"]
    return models.findEFM(model, include, exclude)
    
def emfG():
    include = {"Pfk": models.FORWARD}
    exclude = []
    return models.findEFM(model, include, exclude)
    
def ru5p_semi():
    include = ["Ru5P"]
    testPsemi(include)

if __name__ == '__main__':
    fluxes = {}
    fluxes["A"] = [emfA()]
    fluxes["B"] = [emfB()]
    fluxes["C"] = [emfC()]
    fluxes["D"] = [emfD()]
    fluxes["E"] = [emfE()]
    fluxes["F"] = [emfF()]
    fluxes["G"] = [emfG()]
    display.displayAll("toyModelMap.json", fluxes, "EMFs according to the paper")