# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2018 replay file
# Internal Version: 2017_11_07-19.21.41 127140
# Run by VJAX09 on Wed Jun 20 10:59:29 2018
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=226.170822143555, 
    height=161.977783203125)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='test_dcoupling_brick.odb')
#: Model: C:/Users/vjax09/.julia/v0.6/FEMCoupling/test/test_dcoupling_brick/test_dcoupling_brick_abaqus/test_dcoupling_brick.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       3
#: Number of Node Sets:          6
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    DEFORMED, ))
