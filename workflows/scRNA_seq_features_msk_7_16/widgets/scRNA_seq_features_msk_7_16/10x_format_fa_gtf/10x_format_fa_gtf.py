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

class OW10x_format_fa_gtf(OWBwBWidget):
    name = "10x_format_fa_gtf"
    description = "Formats downloaded gtf fa files in the manner of 10x"
    priority = 20
    icon = getIconName(__file__,"perl.png")
    want_main_area = False
    docker_image_name = "biodepot/format-fa-gtf-10x-sc"
    docker_image_tag = "latest"
    inputs = [("inputgtf",str,"handleInputsinputgtf"),("inputfa",str,"handleInputsinputfa"),("trigger",str,"handleInputstrigger"),("outputgtf",str,"handleInputsoutputgtf"),("outputfa",str,"handleInputsoutputfa"),("skip",str,"handleInputsskip"),("overwritefiles",str,"handleInputsoverwritefiles")]
    outputs = [("outputfa",str),("outputgtf",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    inputfa=pset(None)
    outputfa=pset(None)
    inputgtf=pset(None)
    outputgtf=pset(None)
    overwritefiles=pset(False)
    inputfeaturegtf=pset(None)
    inputfeaturefa=pset(None)
    skip=pset(False)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"10x_format_fa_gtf")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsinputgtf(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("inputgtf", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsinputfa(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("inputfa", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputstrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsoutputgtf(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("outputgtf", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsoutputfa(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("outputfa", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsskip(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("skip", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsoverwritefiles(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("overwritefiles", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputfa"):
            outputValue=getattr(self,"outputfa")
        self.send("outputfa", outputValue)
        outputValue=None
        if hasattr(self,"outputgtf"):
            outputValue=getattr(self,"outputgtf")
        self.send("outputgtf", outputValue)
