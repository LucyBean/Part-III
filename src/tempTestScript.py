'''
Created on Nov 15, 2016

@author: Lucy
'''

import cobra
import models
import escher

externalMetabolites = ["NADP", "NADPH", "CO2", "G6P", "R5Pex", "ATP", "ADP", "Pyr", "NAD", "NADH"]
possibleTerminals = ["G6P", "R5Pex", "Pyr"]
startID = "G6P"
    
cobraModel = cobra.io.load_json_model("toyModel.json")

products = models.findProducts(cobraModel, startID, possibleTerminals, externalMetabolites)

# Make the main visualiser page
startMetaboliteName = "G6P"

with open("../visualisation/working/visualise.html", "wb") as f:
    
    f.write("""<html>
<head>
<style>
iframe {
    width:100%;
    height:100%
}
table {
    height: 100%;
    width: 100%
}
div#visualiserPanel {
    height: 100%
}
/* Set to 1px to make it as small as possible*/
tr#header {
    height: 1px
}
</style>
</head>
<body>
<table>
    <tr id="header">
        <td>
            <h2>Starting metabolite</h2>
            """ + startID + """
            <h2>Possible products</h2>\n""")
    
    index = 1
    for p in products:
        nameIndex = 1
        for flux in products[p]:
            f.write("""\t\t\t<button onclick="showEFM(""" + str(index) + 
                    """)">""" + p + " (" + str(nameIndex) + """)</button>\n""")
            nameIndex += 1
            index += 1
            
    f.write("""\t\t</td>
    </tr>
    <tr>
        <td>
            <div id="visualiserPanel">
            Click on a button to visualise the EFM.
            </div>
        </td>

    </tr>
</table>
<script>
var showEFM = function (num) {
    var displayDiv = document.getElementById("visualiserPanel")
    displayDiv.innerHTML = "<iframe src=\\"efm" + num + ".html\\"></iframe>"
    
}
</script>
</body>
</html>""")

# Output solution HTMLs for iframes
index = 1
for p in products:
    for flux in products[p]:
        b = escher.Builder(map_json="toyModelMap.json", reaction_data=flux,
                           # color and size according to the absolute value
                       reaction_styles=['color', 'size', 'abs', 'text'],
                       # change the default colors
                       reaction_scale=[{'type': 'min', 'color': '#00cc00', 'size': 4},
                                       {'type': 'mean', 'color': '#0000dd', 'size': 20},
                                       {'type': 'max', 'color': '#ff0000', 'size': 40}])
        b.save_html(filepath="../visualisation/working/efm" + str(index) + ".html",
                    overwrite=True, never_ask_before_quit=True, scroll_behavior="zoom")
        index += 1
        

        
