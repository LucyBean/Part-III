'''
Created on Nov 26, 2016

@author: Lucy
'''
import cobra.test
from src import models, display

if __name__ == '__main__':
    startID = "lac__D_e"
    cobraModel = cobra.test.create_test_model("textbook")
    products = models.findProducts(cobraModel, startID)
    title = "Possible products for starting metabolite " + startID + " in textbook model"
    display.displayAll("toyModelMap.json", products, title)