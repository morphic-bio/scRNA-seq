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

class OWGather_features(OWBwBWidget):
    name = "Gather_features"
    description = "Minimum Python3 container with pip"
    priority = 10
    icon = getIconName(__file__,"python3.png")
    want_main_area = False
    docker_image_name = "biodepot/gather_features"
    docker_image_tag = "latest"
    inputs = [("inputFile",str,"handleInputsinputFile"),("Trigger",str,"handleInputsTrigger"),("aligndir",str,"handleInputsaligndir"),("Trigger2",str,"handleInputsTrigger2")]
    outputs = [("OutputDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    features_dirs=pset([])
    features_gex_dirs=pset(None)
    aligndir=pset(None)
    cellbender_file=pset(None)
    cb_original_layer=pset(None)
    cb_final_layer=pset(None)
    input_counts_name=pset("unfiltered_counts.h5ad")
    output_features_name=pset("merged_features.h5ad")
    output_counts_name=pset("merged_counts.h5ad")
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"Gather_features")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsinputFile(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("inputFile", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsTrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("Trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsaligndir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("aligndir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsTrigger2(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("Trigger2", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"OutputDir"):
            outputValue=getattr(self,"OutputDir")
        self.send("OutputDir", outputValue)
