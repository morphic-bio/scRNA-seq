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

class OWfilterY(OWBwBWidget):
    name = "filterY"
    description = "Remove Y chromosome reads"
    priority = 10
    icon = getIconName(__file__,"filtered_y_chromosome_resized_square.png")
    want_main_area = False
    docker_image_name = "biodepot/yfilter"
    docker_image_tag = "latest"
    inputs = [("trigger",str,"handleInputstrigger")]
    outputs = [("outputDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    aligneddir=pset(None)
    sequenceDir=pset(None)
    outputDir=pset(None)
    overwritefiles=pset(False)
    nThreads=pset(4)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"filterY")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputstrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputDir"):
            outputValue=getattr(self,"outputDir")
        self.send("outputDir", outputValue)
