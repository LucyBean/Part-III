'''
Created on Nov 11, 2016

@author: Lucy
'''

import cobra.test
import os
from os.path import join

data_dir = cobra.test.data_directory

print "Mini test files:"
print ", ".join(i for i in os.listdir(data_dir) if i.startswith("mini"))

textbook_model = cobra.test.create_test_model("textbook")
ecoli_model = cobra.test.create_test_model("ecoli")
salmonella_model = cobra.test.create_test_model("salmonella")
mini = cobra.io.read_sbml_model(join(data_dir, "mini_fbc2.xml"))
iJO1366 = cobra.io.read_sbml_model("MODEL1108160000.xml")
 
