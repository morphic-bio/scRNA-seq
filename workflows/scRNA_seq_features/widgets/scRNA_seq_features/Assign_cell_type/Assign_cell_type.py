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

class OWAssign_cell_type(OWBwBWidget):
    name = "Assign_cell_type"
    description = "Enter and output a file"
    priority = 10
    icon = getIconName(__file__,"Transformer.png")
    want_main_area = False
    docker_image_name = "biodepot/alpine-bash"
    docker_image_tag = "3.7"
    inputs = [("File",str,"handleInputsFile")]
    outputs = [("File",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    filename=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"Assign_cell_type")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsFile(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("File", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"File"):
            outputValue=getattr(self,"File")
        self.send("File", outputValue)
