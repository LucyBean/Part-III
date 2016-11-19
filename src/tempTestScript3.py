'''
Created on Nov 15, 2016

@author: Lucy
'''
reactants = []
products = []

# Extract reactants and products from reactions
for f in flux:
    reaction = reactions.get_by_id(f)
    rf = reaction.reactants
    pf = reaction.products
    if flux[f] > 0:
        reactants += [r.id for r in rf]
        products += [p.id for p in pf]
    else:
        reactants += [p.id for p in pf]
        products += [r.id for r in rf]
        
# extract the metabolites that are only produced or consumed
onlyProduced = [m for m in reactants if m not in products]
onlyConsumed = [m for m in products if m not in reactants]