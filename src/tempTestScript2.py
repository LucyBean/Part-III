flux = models.findEFM(cobraModel, reactionsToInclude, reactionsToExclude, externalMetabolites)

if flux is None:
    print "Flux is none."

possibleTerminals = ["G6P", "R5Pex", "Pyr"]

(terminalReactants, terminalProducts) = models.findTerminalReactantsAndProducts(possibleTerminals, flux, reactions)

# Take the first product in terminalProducts and knock out that reaction
if terminalProducts != []:
    terminalProduct = metabolites.get_by_id(terminalProducts[0])
    rs = list(terminalProduct.reactions)
    # Take the first reaction in this list
    knockOut = [r for r in rs if r.id in flux][0]
    reactionsToExclude.append(knockOut.id)
# The EMF produced is a loop
else:
    terminalProduct = startingMetabolite
    rs = list(terminalProduct.reactions)
    # Find all the reactions that the starting metabolite is involved in
    mightKnockOut = [r for r in rs if r.id in flux]
    # Find all of the ones that produce the starting metabolite
    knockOut = []
    for r in mightKnockOut:
        if flux[r.id] > 0 and startingMetabolite in r.products:
            knockOut.append(r)
        elif flux[r.id] < 0 and startingMetabolite in r.reactants:
            knockOut.append(r)
    reactionsToExclude.append(knockOut[0].id)