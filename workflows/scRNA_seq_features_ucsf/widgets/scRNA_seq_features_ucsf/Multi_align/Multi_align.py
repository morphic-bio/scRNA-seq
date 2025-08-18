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

class OWMulti_align(OWBwBWidget):
    name = "Multi_align"
    description = "Construct indices for STAR aligner "
    priority = 11
    icon = getIconName(__file__,"staralign.png")
    want_main_area = False
    docker_image_name = "biodepot/multialign"
    docker_image_tag = "latest"
    inputs = [("Trigger",str,"handleInputsTrigger"),("fastqdirs",str,"handleInputsfastqdirs"),("outputdir",str,"handleInputsoutputdir"),("indicesdir",str,"handleInputsindicesdir"),("useStar",str,"handleInputsuseStar"),("useSalmon",str,"handleInputsuseSalmon"),("useKallisto",str,"handleInputsuseKallisto"),("usePiscem",str,"handleInputsusePiscem"),("useSplicei",str,"handleInputsuseSplicei"),("useSpliceu",str,"handleInputsuseSpliceu"),("skip",str,"handleInputsskip")]
    outputs = [("outputdir",str),("sampledirs",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    fastqdirs=pset([])
    outputdir=pset(None)
    ntreads=pset(8)
    indicesdir=pset(None)
    useStar=pset(False)
    useKallisto=pset(False)
    usePiscem=pset(False)
    useSalmon=pset(False)
    useSplicei=pset(False)
    useSpliceu=pset(False)
    starcmd=pset(None)
    kallistocmd=pset(None)
    salmoncmd=pset(None)
    piscemcmd=pset(None)
    overwrite=pset(False)
    whitelist=pset(None)
    R1pattern=pset("_R1_")
    R2pattern=pset("_R2_")
    gzipped_fastq=pset(False)
    fastq_suffix=pset("fastq")
    CBLEN=pset(16)
    UMILEN=pset(12)
    organize_samples=pset(False)
    concatenate_files=pset(False)
    sampledirs=pset([])
    skip=pset(False)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"Multi_align")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsTrigger(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("Trigger", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsfastqdirs(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("fastqdirs", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsoutputdir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("outputdir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsindicesdir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("indicesdir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsuseStar(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("useStar", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsuseSalmon(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("useSalmon", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsuseKallisto(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("useKallisto", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsusePiscem(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("usePiscem", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsuseSplicei(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("useSplicei", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsuseSpliceu(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("useSpliceu", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsskip(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("skip", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputdir"):
            outputValue=getattr(self,"outputdir")
        self.send("outputdir", outputValue)
        outputValue=None
        if hasattr(self,"sampledirs"):
            outputValue=getattr(self,"sampledirs")
        self.send("sampledirs", outputValue)
