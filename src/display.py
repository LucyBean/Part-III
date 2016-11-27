'''
Created on Nov 24, 2016

@author: Lucy
'''

import escher
import os
import datetime

def displayEFM(map_name=None, map_json=None, reaction_data=[]):
    """Displays an EFM on the given map using the given reaction data."""
    b = escher.Builder(map_name=map_name, map_json=map_json,
                       reaction_data=reaction_data,
                       # color and size according to the absolute value
                   reaction_styles=['color', 'size', 'abs', 'text'],
                   # change the default colors
                   reaction_scale=[{'type': 'min', 'color': '#009933', 'size': 4},
                                   {'type': 'mean', 'color': '#0000dd', 'size': 20},
                                   {'type': 'max', 'color': '#ff0000', 'size': 40}])
    dirID = datetime.datetime.now().strftime("%Y-%m-%d %H%M%S")
    os.makedirs("visualisation/" + dirID)
    filepath = "visualisation/" + dirID + "/efm.html"
    b.save_html(filepath=filepath, overwrite=True, never_ask_before_quit=True, scroll_behavior="zoom")
    path = os.path.abspath(filepath)
    os.startfile(path)
    
def displayPsemi(map_name=None, map_json=None, metabolite_data=[]):
    """Displays a P-semiflow on the given map using the given metabolite data."""
    b = escher.Builder(map_name=map_name, map_json=map_json,
                       metabolite_data=metabolite_data,
                       # color and size according to the absolute value
                   reaction_styles=['color', 'size', 'abs', 'text'],
                   # change the default colors
                   reaction_scale=[{'type': 'min', 'color': '#009933', 'size': 4},
                                   {'type': 'mean', 'color': '#0000dd', 'size': 20},
                                   {'type': 'max', 'color': '#ff0000', 'size': 40}])
    dirID = datetime.datetime.now().strftime("%Y-%m-%d %H%M%S")
    os.makedirs("visualisation/" + dirID)
    filepath = "visualisation/" + dirID + "/efm.html"
    b.save_html(filepath=filepath, overwrite=True, never_ask_before_quit=True, scroll_behavior="zoom")
    path = os.path.abspath(filepath)
    os.startfile(path)
    
def displayAll(map_json, toDisplay, title=""):
    """Displays all of the EFMs in toDisplay for the given cobraModel, map_json
    
    toDisplay: A dict of the form {name: list of fluxes}. Each flux will be have a button
                on the generated web page which can be used to display the corresponding
                EFM. This button will be have a label "name (#)" where # is a number.
    """
    dirID = datetime.datetime.now().strftime("%Y-%m-%d %H%M%S") + " " + title
    os.makedirs("visualisation/" + dirID)
    with open("visualisation/" + dirID + "/visualise.html", "wb") as f:
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
        <td><h1>""" + title + "</h1>")
        
        index = 1
        for p in toDisplay:
            nameIndex = 1
            for flux in toDisplay[p]:
                # Put a button on the main page
                f.write("""\t\t\t<button onclick="showEFM(""" + str(index) + 
                        """)">""" + p)
                if len(toDisplay[p]) != 1:
                    f.write(" (" + str(nameIndex) + ")")
                    nameIndex += 1
                f.write("</button>\n")
                
                # Output this flux as an HTML file
                b = escher.Builder(map_json=map_json, reaction_data=flux,
                                   # color and size according to the absolute value
                               reaction_styles=['color', 'size', 'abs', 'text'],
                               # change the default colors
                               reaction_scale=[{'type': 'min', 'color': '#009933', 'size': 4},
                                               {'type': 'mean', 'color': '#0000dd', 'size': 20},
                                               {'type': 'max', 'color': '#ff0000', 'size': 40}])
                b.save_html(filepath="visualisation/" + dirID + "/efm" + str(index) + ".html",
                            overwrite=True, never_ask_before_quit=True, scroll_behavior="zoom")
                
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
    
    # Display the file
    path = os.path.abspath(f.name)
    os.startfile(path)  
    
