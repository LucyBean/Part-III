'''
Created on Nov 14, 2016

@author: Lucy
'''

import models
import cobra

defaultIgnore = ["NADP", "NADPH", "CO2", "G6P", "R5Pex", "ATP", "ADP", "Pyr", "NAD", "NADH"]

def test(include, exclude = [], ignore=defaultIgnore):
    model = cobra.io.load_json_model("toyModel.json")
    
    fluxes = models.process(model, include, exclude, ignore)
    if fluxes is not None:
        models.display(map_json="toyModelMap.json", reaction_data=fluxes)
        
def emfA():
    include = {"Pgi": models.FORWARD, "Pyk": models.FORWARD}
    exclude = ["Prs_DeoB", "Zwf"]
    test(include, exclude)
    
def emfB():
    include = {"Zwf": models.FORWARD}
    exclude = ["Prs_DeoB", "Pfk", "Fbp"]
    test(include, exclude)
    
def emfC():
    include = {"Zwf": models.FORWARD}
    exclude = ["Prs_DeoB", "Pgi"]
    test(include, exclude)
    
def emfD():
    include = {"Zwf": models.FORWARD}
    exclude = []
    test(include, exclude)
    
def emfE():
    include = {"Pgi": models.FORWARD}
    exclude = []
    test(include, exclude)
    
def emfF():
    include = {"Zwf": models.FORWARD}
    exclude = ["Prs_DeoB"]
    test(include, exclude)
    
def emfG():
    include = {"Pfk": models.FORWARD}
    exclude = []
    test(include, exclude)

if __name__ == '__main__':
    emfB()