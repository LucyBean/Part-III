'''
Created on Nov 14, 2016

@author: Lucy
'''

from cobra import Model, Reaction, Metabolite
from cobra.io import save_json_model

def reversibleReaction(name):
    return Reaction(name, lower_bound=-1000)

if __name__ == '__main__':
    model = Model("Monosaccharide metabolism")
    
    # Make metabolites
    G6P = Metabolite("G6P")
    F6P = Metabolite("F6P")
    FP2 = Metabolite("FP2")
    GAP = Metabolite("GAP")
    BPG13 = Metabolite("1.3BPG")
    PG3 = Metabolite("3PG")
    PG2 = Metabolite("PG2")
    PEP = Metabolite("PEP")
    PYR = Metabolite("Pyr")
    Ery4P = Metabolite("Ery4P")
    Sed7P = Metabolite("Sed7P")
    R5P = Metabolite("R5P")
    R5Pex = Metabolite("R5Pex")
    Xyl5P = Metabolite("Xyl5P")
    Ru5P = Metabolite("Ru5P")
    PG6 = Metabolite("6PG")
    GO6P = Metabolite("GO6P")
    DHAP = Metabolite("DHAP")
    NADP = Metabolite("NADP")
    NADPH = Metabolite("NADPH")
    ADP = Metabolite("ADP")
    ATP = Metabolite("ATP")
    NAD = Metabolite("NAD")
    NADH = Metabolite("NADH")
    CO2 = Metabolite("CO2")
    
    # Make reactions
    Pgi = reversibleReaction("Pgi")
    Zwf = reversibleReaction("Zwf")
    Pgl = reversibleReaction("Pgl")
    Gnd = reversibleReaction("Gnd")
    Rpi = reversibleReaction("Rpi")
    Rpe = reversibleReaction("Rpe")
    Prs_DeoB = Reaction("Prs_DeoB")
    TktI = reversibleReaction("TktI")
    Tal = reversibleReaction("Tal")
    TktII = reversibleReaction("TktII")
    Fbp = Reaction("Fbp")
    Pfk = Reaction("Pfk")
    Fba = reversibleReaction("Fba")
    TpiA = reversibleReaction("TpiA")
    Gap = reversibleReaction("Gap")
    Pgk = reversibleReaction("Pgk")
    Gpm = reversibleReaction("Gpm")
    Eno = reversibleReaction("Eno")
    Pyk = Reaction("Pyk")
    
    # Add metabolites to reactions
    Pgi.add_metabolites({G6P:-1, F6P: 1})
    Zwf.add_metabolites({G6P:-1, NADP:-1, NADPH: 1, GO6P: 1})
    Pgl.add_metabolites({GO6P:-1, PG6: 1})
    Gnd.add_metabolites({PG6:-1, NADP:-1, NADPH: 1, CO2: 1, Ru5P: 1})
    Rpi.add_metabolites({Ru5P:-1, R5P: 1})
    Rpe.add_metabolites({Ru5P:-1, Xyl5P: 1})
    Prs_DeoB.add_metabolites({R5P:-1, R5Pex: 1})
    TktI.add_metabolites({Xyl5P:-1, R5P:-1, Sed7P: 1, GAP: 1})
    Tal.add_metabolites({Sed7P:-1, GAP:-1, Ery4P: 1, F6P: 1})
    TktII.add_metabolites({Xyl5P:-1, Ery4P:-1, F6P: 1, GAP: 1})
    Fbp.add_metabolites({FP2:-1, F6P: 1})
    Pfk.add_metabolites({F6P:-1, ATP:-1, FP2: 1, ADP: 1})
    Fba.add_metabolites({FP2:-1, GAP: 1, DHAP: 1})
    TpiA.add_metabolites({DHAP:-1, GAP: 1})
    Gap.add_metabolites({GAP:-1, NAD:-1, NADH: 1, BPG13: 1})
    Pgk.add_metabolites({BPG13:-1, ADP:-1, ATP: 1, PG3: 1})
    Gpm.add_metabolites({PG3:-1, PG2: 1})
    Eno.add_metabolites({PG2:-1, PEP: 1})
    Pyk.add_metabolites({PEP:-1, ADP:-1, ATP: 1, PYR: 1})
    
    model.add_reactions([Pgi, Zwf, Pgl, Gnd, Rpi, Rpe, Prs_DeoB, TktI,
                         Tal, TktII, Fbp, Pfk, Fba, TpiA, Gap, Pgk,
                         Gpm, Eno, Pyk])
    
    save_json_model(model, "toyModel.json")
