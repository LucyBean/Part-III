# hack hack
if "products" not in locals():
    execfile("tempTestScript.py")

startMetaboliteName = "G6P"

with open("../visualisation/test/test.html", "wb") as f:
    
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
td {
    background-color:#7777FF;
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
            """ + startMetaboliteName + """
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
    