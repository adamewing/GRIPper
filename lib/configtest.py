#!/bin/env python

"""
 Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)

 Released under the MIT license, see LICENSE.txt

 configtest.py checks required configuration options and files
"""

import ConfigParser
import sys
import os

def checkfile(file):
    try:
        open(file)
        return True
    except:
        return False

def check(cfgfile):

    if not checkfile(cfgfile):
        raise ValueError("Could not open config file: " + cfgfile)

    valoptions  = ('readLength',
                   'insertSize',
                   'minMapQ',
                   'minPeakSize',
                   'maxReadLen',
                   'zeroChar',
                   'minClipQual',
                   'configFileName')

    diroptions  = ('tabixDir',
                   'annotDir')

    fileoptions = ('pgTabixFile',
                   'mrnaFastaFile')

    config = ConfigParser.ConfigParser()
    config.read(cfgfile)

    for option in valoptions:
        if not config.has_option('discord', option):
            raise ValueError("option " + option + " not found.")

    for option in fileoptions:
        if not config.has_option('discord', option):
            raise ValueError("option " + option + " not found.")
        filename = config.get('discord', option)
        if not checkfile(filename):
            raise ValueError("file not found: " + filename + " for option: " + option) 

    for option in diroptions:
        dir = config.get('discord', option)
        if not os.path.exists(dir):
            raise ValueError("directory not found: " + dir)
 
if __name__ == '__main__':
    if len(sys.argv) == 2:
        check(sys.argv[1])
