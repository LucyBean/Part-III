'''
Created on Nov 12, 2016

@author: Lucy
'''

from cobra import Model, Reaction, Metabolite
    
cobra_model = Model("example_cobra_model")


# ## Make the reaction
reaction = Reaction("3OAS140")
reaction.name = "3 oxyacyl acyl carrier protein synthase n C140"
reaction.subsystem = "Cell Envelope Biosynthesis"
reaction.lower_bound = 0  # default
reaction.upper_bound = 1000  # default
reaction.objective_coefficient = 0  # default


# ## Make the metabolites
ACP_c = Metabolite(
    'ACP_c',
    formula='C11H21N2O7PRS',
    name='acyl-carrier-protein',
    compartment='c')
omrsACP_c = Metabolite(
    '3omrsACP_c',
    formula='C25H45N2O9PRS',
    name='3-Oxotetradecanoyl-acyl-carrier-protein',
    compartment='c')
co2_c = Metabolite(
    'co2_c',
    formula='CO2',
    name='CO2',
    compartment='c')
malACP_c = Metabolite(
    'malACP_c',
    formula='C14H22N2O10PRS',
    name='Malonyl-acyl-carrier-protein',
    compartment='c')
h_c = Metabolite(
    'h_c',
    formula='H',
    name='H',
    compartment='c')
ddcaACP_c = Metabolite(
    'ddcaACP_c',
    formula='C23H43N2O8PRS',
    name='Dodecanoyl-ACP-n-C120ACP',
    compartment='c')


# ## Add metabolites to the reaction
reaction.add_metabolites({malACP_c:-1.0,
                          h_c:-1.0,
                          ddcaACP_c:-1.0,
                          co2_c: 1.0,
                          ACP_c: 1.0,
                          omrsACP_c: 1.0})


# ## Genes
reaction.gene_reaction_rule = "(STM2378 or STM1197)"


# ## Add reaction to the model
cobra_model.add_reaction(reaction)
print "Reactions"
for x in cobra_model.reactions:
    print "%s : %s" % (x.id, x.reaction)
print "\n"
print "Metabolites"
for x in cobra_model.metabolites:
    print "%s : %s" % (x.id, x.formula)
print "\n"
print "Genes"
for x in cobra_model.genes:
    print x 
