# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2020 replay file
# Internal Version: 2019_09_13-13.49.31 163176
# Run by vtac on Wed Dec 15 20:57:29 2021
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=368.979156494141, 
    height=176.308944702148)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
o1 = session.openOdb(name='/home/vtac/NODE_abaqus/UMAT/uni.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/vtac/NODE_abaqus/UMAT/uni.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       6
#: Number of Node Sets:          6
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
session.printOptions.setValues(vpDecorations=OFF)
session.printToFile(fileName='/home/vtac/NODE_abaqus/images/UMA_uni1.png', 
    format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.printToFile(fileName='/home/vtac/NODE_abaqus/images/UMA_uni2.png', 
    format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))
o1 = session.openOdb(name='/home/vtac/NODE_abaqus/UMAT/tor.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/vtac/NODE_abaqus/UMAT/tor.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             2
#: Number of Element Sets:       6
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.printToFile(fileName='/home/vtac/NODE_abaqus/images/UMA_tor1.png', 
    format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.printToFile(fileName='/home/vtac/NODE_abaqus/images/UMA_tor2.png', 
    format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))
o1 = session.openOdb(name='/home/vtac/NODE_abaqus/UMAT/she.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/vtac/NODE_abaqus/UMAT/she.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       8
#: Number of Node Sets:          8
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.printToFile(fileName='/home/vtac/NODE_abaqus/images/UMA_she1.png', 
    format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.printToFile(fileName='/home/vtac/NODE_abaqus/images/UMA_she2.png', 
    format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))
o1 = session.openOdb(name='/home/vtac/NODE_abaqus/UMAT/exp.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/vtac/NODE_abaqus/UMAT/exp.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     3
#: Number of Meshes:             4
#: Number of Element Sets:       3
#: Number of Node Sets:          5
#: Number of Steps:              3
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.printToFile(fileName='/home/vtac/NODE_abaqus/images/UMA_exp1.png', 
    format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))
