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
# use sets to prevent extra work when handling duplicates
onlyProduced = set([m for m in reactants if m not in products])
onlyConsumed = set([m for m in products if m not in reactants])

# find the reaction(s?) that consumed them
for m in products:
    metabolite = metabolites.get_by_id(m)
    # Find the reaction(s?) that include this metabolite
    myReactions = [r.id for r in list(metabolite.reactions) if r.id in flux]
    # Find those reactions by which this metabolite is consumed
    # This will be all forward reactions for which this metabolite is a reactant
    # Or all revere reactions for which this metabolite is a product
    myConsumers = []
    for mr in myReactions:
        f = flux[mr]
        mrs = reactions.get_by_id(mr)
        # If this is included as a forward reaction
        if f > 0:
            # Find the IDs of all reactants
            rs = [r.id for r in mrs.reactants]
            # And include it if this metabolite is included
            if m in rs:
                myConsumers.append(mr)
        # If this is included as a reverse reaction
        elif f < 0:
            # Find the IDs of all products
            rs = [r.id for r in mrs.products]
            # And include it if this metabolite is included
            if m in rs:
                myConsumers.append(mr)
        
            
    print metabolite, myConsumers