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

class OWfeature_finder(OWBwBWidget):
    name = "feature_finder"
    description = "Assigns feature barcodes to cells"
    priority = 13
    icon = getIconName(__file__,"Barcode_font_awesome.svg.png")
    want_main_area = False
    docker_image_name = "biodepot/process_features"
    docker_image_tag = "latest"
    inputs = [("outputdir",str,"handleInputsoutputdir"),("trigger",str,"handleInputstrigger"),("features_file",str,"handleInputsfeatures_file"),("whitelist",str,"handleInputswhitelist")]
    outputs = [("outputdir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    imputfiles=pset([])
    inputdirs=pset([])
    keep_existing=pset(False)
    whitelist=pset(None)
    barcode_length=pset(16)
    features_file=pset(None)
    maxHammingDistance=pset(5)
    stringency=pset(1)
    min_counts=pset(0)
    outputdir=pset("")
    feature_constant_offset=pset(26)
    nprocesses=pset(None)
    searchThreads=pset(4)
    consumerThreads=pset(None)
    maxThreads=pset(None)
    parallelbyFile=pset(False)
    umi_length=pset(12)
    debugmode=pset(False)
    process_in_order=pset(False)
    maximum_mismatch=pset(3)
    barcode_n=pset(1)
    feature_n=pset(3)
    barcode_constant_offset=pset(None)
    readBufferLines=pset(2048)
    as_named=pset(False)
    reverse_complement=pset(False)
    barcodePattern=pset("_R1_")
    forwardPattern=pset("_R2_")
    reversePattern=pset("_R3_")
    averageReadLength=pset(None)
    minPosterior=pset(0.975)
    skip=pset(False)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"feature_finder")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsoutputdir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("outputdir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputstrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsfeatures_file(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("features_file", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputswhitelist(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("whitelist", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputdir"):
            outputValue=getattr(self,"outputdir")
        self.send("outputdir", outputValue)
