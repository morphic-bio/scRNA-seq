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

class OWGather_scRNA_seq(OWBwBWidget):
    name = "Gather_scRNA_seq"
    description = "alpine bash with wget curl gzip bzip2"
    priority = 1
    icon = getIconName(__file__,"bash.png")
    want_main_area = False
    docker_image_name = "biodepot/gather_rna_seq"
    docker_image_tag = "latest"
    inputs = [("counts_dir",str,"handleInputscounts_dir"),("Trigger1",str,"handleInputsTrigger1")]
    outputs = [("output_dir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    counts_dir=pset([])
    counts_patterns=pset([])
    features_dirs=pset([])
    aligns_dirs=pset([])
    features_patterns=pset(None)
    aligns_patterns=pset([])
    output_dir=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"Gather_scRNA_seq")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputscounts_dir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("counts_dir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsTrigger1(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("Trigger1", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"output_dir"):
            outputValue=getattr(self,"output_dir")
        self.send("output_dir", outputValue)
