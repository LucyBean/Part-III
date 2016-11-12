'''
Created on Nov 12, 2016

@author: Lucy
'''
import escher
import cobra

model = cobra.io.json.from_json(escher.plots.model_json_for_name('e_coli_core'))
solution = model.optimize()
b = escher.Builder(map_name='e_coli_core.Core metabolism',
                   reaction_data=solution.x_dict,
                   # color and size according to the absolute value
                   reaction_styles=['color', 'size', 'abs', 'text'],
                   # change the default colors
                   reaction_scale=[{'type': 'min', 'color': '#cccccc', 'size': 4},
                                   {'type': 'mean', 'color': '#0000dd', 'size': 20},
                                   {'type': 'max', 'color': '#ff0000', 'size': 40}],
                   # only show the primary metabolites
                   hide_secondary_metabolites=True)
b.display_in_browser(scroll_behavior="zoom")