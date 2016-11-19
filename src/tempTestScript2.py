'''
Created on Nov 15, 2016

@author: Lucy
'''



# extract the metabolites that are only produced or consumed
# use sets to prevent extra work when handling duplicates
onlyProduced = set([m for m in products if m not in reactants])
onlyConsumed = set([m for m in reactants if m not in products])

terminalReactants = []
# find the reaction(s?) that consumed them
for m in onlyConsumed:
    metabolite = metabolites.get_by_id(m)
    # Find the reaction(s?) that include this metabolite
    myReactions = [r.id for r in list(metabolite.reactions) if r.id in flux]
    # Find those reactions by which this metabolite is consumed
    # This will be all forward reactions for which this metabolite is a reactant
    # Or all reverse reactions for which this metabolite is a product
    print m, myReactions
    
    for mrid in myReactions:
        f = flux[mrid]
        mrr = reactions.get_by_id(mrid)
        # Need to find whether this reactant is the ONLY reactant, in which case
        #  it is definitely a terminal reactant
        # Get the list of reactants for the reaction. This depends on the direction
        #  of the reaction, which is determined from the flux vector.
        rs = []
        if f > 0:
            # For a forward reaction, reactants are as expected
            rs = [r.id for r in mrr.reactants]
        elif f < 0:
            # For a backwards reaction, the reactants are the reaction's products
            rs = [r.id for r in mrr.products]
        
        # If this is the only reactant then it must be a terminal reactant
        if len(rs) == 1:
            terminalReactants += [m]
        
print terminalReactants
        
    