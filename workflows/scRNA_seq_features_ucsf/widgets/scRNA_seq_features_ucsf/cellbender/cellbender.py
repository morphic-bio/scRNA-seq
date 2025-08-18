import os
import glob
import sys
import functools
import jsonpickle
from collections import OrderedDict
from Orange.widgets import widget, gui, settings
import Orange.data
from Orange.data.io import FileFormat
from DockerClient import DockerClient
from BwBase import OWBwBWidget, ConnectionDict, BwbGuiElements, getIconName, getJsonName
from PyQt5 import QtWidgets, QtGui

class OWcellbender(OWBwBWidget):
    name = "cellbender"
    description = "Gathers counts from star/salmon workflow"
    priority = 13
    icon = getIconName(__file__,"logo_250_185.png")
    want_main_area = False
    docker_image_name = "biodepot/cellbender"
    docker_image_tag = "0.3.2"
    inputs = [("countsDir",str,"handleInputscountsDir"),("alignsDir",str,"handleInputsalignsDir"),("tablesDir",str,"handleInputstablesDir"),("trigger",str,"handleInputstrigger")]
    outputs = [("tablesDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    alignsDir=pset([])
    overwrite_cellbender=pset(False)
    input_pattern=pset("filter_counts.h5ad")
    output_pattern=pset("final_counts.h5ad")
    nThreads=pset(1)
    cb_subdir=pset("cellbender")
    cb_counts_file=pset("cellbender_counts.h5")
    overwrite_layer=pset(False)
    usecpu=pset(False)
    layername=pset("denoised")
    cpu_cores=pset(1)
    additional_flags=pset(None)
    skip=pset(False)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"cellbender")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputscountsDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("countsDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsalignsDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("alignsDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputstablesDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("tablesDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputstrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"tablesDir"):
            outputValue=getattr(self,"tablesDir")
        self.send("tablesDir", outputValue)
