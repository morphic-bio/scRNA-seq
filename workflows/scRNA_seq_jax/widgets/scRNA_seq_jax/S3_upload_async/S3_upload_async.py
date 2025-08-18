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

class OWS3_upload_async(OWBwBWidget):
    name = "S3_upload_async"
    description = "Enter and output a file"
    priority = 10
    icon = getIconName(__file__,"cloud_upload.png")
    want_main_area = False
    docker_image_name = "biodepot/s3upload_async"
    docker_image_tag = "latest"
    inputs = [("Trigger",str,"handleInputsTrigger"),("credentials_dir",str,"handleInputscredentials_dir"),("uploadDir",str,"handleInputsuploadDir"),("bucket",str,"handleInputsbucket"),("s3Dir",str,"handleInputss3Dir")]
    outputs = [("credentials_dir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    credentials_dir=pset("/data/.aws")
    uploadDir=pset([])
    bucket=pset(None)
    s3Dir=pset(None)
    profile=pset(None)
    doneDir=pset(None)
    sleeptime=pset(60)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"S3_upload_async")) as f:
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
    def handleInputscredentials_dir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("credentials_dir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsuploadDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("uploadDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsbucket(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("bucket", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputss3Dir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("s3Dir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"credentials_dir"):
            outputValue=getattr(self,"credentials_dir")
        self.send("credentials_dir", outputValue)
