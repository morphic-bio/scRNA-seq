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

class OWFilter_counts(OWBwBWidget):
    name = "Filter_counts"
    description = "Gathers counts from star/salmon workflow"
    priority = 13
    icon = getIconName(__file__,"startodeseq2.png")
    want_main_area = False
    docker_image_name = "biodepot/scrna-matrices"
    docker_image_tag = "latest"
    inputs = [("alignsDir",str,"handleInputsalignsDir"),("trigger",str,"handleInputstrigger")]
    outputs = [("tablesDir",str),("alignsDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    alignsDir=pset([])
    overwrite=pset(False)
    mingenes=pset(200)
    maxgenes=pset(2500)
    mtpercent=pset(5.0)
    nthreads=pset(1)
    skip=pset(False)
    n_mad=pset(3)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"Filter_counts")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsalignsDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("alignsDir", value, args[0][0], test=args[0][3])
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
        outputValue=None
        if hasattr(self,"alignsDir"):
            outputValue=getattr(self,"alignsDir")
        self.send("alignsDir", outputValue)
