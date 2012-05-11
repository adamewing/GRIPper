#!/bin/env python

"""
 Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)

 Released under the MIT license, see LICENSE.txt
"""

import pysam
import argparse
import pairinfo
import pickreads 
import sys
import os
import math
import logging

class PeakBuilder:
    """facilitates keeping reads from different element classes seperate"""
    def __init__(self,geneName):
        self.peak     = pairinfo.Peak()
        self.peaks    = []
	self.geneName = geneName
        self.lastPeakNum = 0

def checkOutDir(outBaseName,outDirName):
    if not os.path.exists(outDirName + "/" + outBaseName):
        raise IOError('cannot find output directory (' + outDirName + ')')
    if not os.path.exists(outDirName + "/" + outBaseName):
        raise IOError("cannot find sample directory (" + outDirName + "/" + outBaseName + ")")
    if not os.path.exists(outDirName + "/" + outBaseName + '/logs'):
        os.mkdir(outDirName + "/" + outBaseName + '/logs')

def checkPseudo(tabixFile,chr,start,end):
    """
    Function for checking coordinates against pseudogene tabix file
    Returns True if an element overlapped
    """
    if start > end: # can happen in messy invalid peaks
       tmp   = start
       start = end
       end   = tmp
    tabixTupleParse = tabixFile.fetch(reference=chr, start=start, end=end, parser=pysam.asTuple())
    if tabixTupleParse:
        for tabixTuple in tabixTupleParse:
            return True
    return False

def readConfig(configPath,outBaseName,outDirName):
    """reads and validates config.txt from output directory"""
    f = open(configPath, 'r')
    configDict = dict()
    for configline in f:
        (param, value) = configline.rstrip().rsplit("=")
        configDict[param] = value

    # check for minimum necessary parameters and some sanity checks
    if configDict['outBaseName']:
        if configDict['outBaseName'] != outBaseName:
            raise ValueError("baseName mismatch: " + outBaseName + "=/="
                             + configDict['outBaseName'])
    else:
        raise ValueError("outBaseName is not set in " + configPath)
    if not configDict['exonTabixFile']:
        raise ValueError("exonTabixFile is not set in " + configPath)

    configDict['readFileName'] = outDirName + "/" + outBaseName + "/" + outBaseName + ".readpairs.txt"

    return configDict

def longOutput(geneNameDict,outBaseName,outDirName,pgTabix):
    """
    Outputs information for the debug.txt file in the output directory
    which contains information for downstream tasks
    """
    fname = outBaseName + ".debug.txt"
    f = open(outDirName + "/" + outBaseName + "/" + fname, 'w')

    for geneName in geneNameDict.keys():
        peakList = geneNameDict[geneName]
        for peak in peakList.peaks:
            if peak.hasPairs():
                breakpoint = peak.getBreakpoint()
		pgOverlap = checkPseudo(pgTabix,peak.getChr(),breakpoint.minPeakPos(),breakpoint.maxPeakPos())
                primaryGene = peak.primaryMateElt().split('.')[0]

                outstr = ("pos=" + peak.getChr() + ":%d-%d" % (breakpoint.leftpos,breakpoint.rightpos)
                        + " str=" + breakpoint.leftstrand + "," + breakpoint.rightstrand + " valid=%d" % breakpoint.isValid()
                        + " lmd=%d" % breakpoint.leftMeanDist() + " rmd=%d" % breakpoint.rightMeanDist()
                        + " lws=%d" % breakpoint.leftWrongStrand() + " rws=%d" % breakpoint.rightWrongStrand()
                        + " lMws=%d" % breakpoint.leftMateWrongStrand() + " rMws=%d" % breakpoint.rightMateWrongStrand()
                        + " np=%d" % breakpoint.npairs + " intD=%d" % (breakpoint.rightpos-breakpoint.leftpos)
                        + " pME=" + peak.primaryMateElt() + " nEx=%d" % peak.numExons() + " pgO=" + str(pgOverlap)
                        + " eltExt=%d,%d" % (breakpoint.eltMinPos(),breakpoint.eltMaxPos())
                        + " exT=" + ",".join(peak.exonTypes()) + " pGene=" + primaryGene
                        + " lExT=" + breakpoint.leftMajorExonType() + " rExT=" + breakpoint.rightMajorExonType()
                        + " eS=" + breakpoint.eltStrandFromReads() + " sources=" + ','.join(breakpoint.sourceList())
                        + " eInv=" + breakpoint.eltStrandFromStrand() + " index=%d" % peak.getIndex() + "\n")
                f.write(outstr)
    f.close()

def main(args):
    checkOutDir(args.outBaseName,args.outDirName)
    configPath = args.outDirName + "/" + args.outBaseName + "/" + args.configFileName
    pickreads.checkfile(configPath)
    configDict = readConfig(configPath,args.outBaseName,args.outDirName)

    # debugging info
    logfile = args.outDirName + "/" + args.outBaseName + "/logs/%d" % os.getpid() + "." + args.outBaseName +".peakparser.log"
    logging.basicConfig(format='%(asctime)s %(message)s',filename=logfile,level=logging.DEBUG)

    # log parameters
    logging.info("\noutBaseName=%s\nbamFileName=%s\nexonTabixFile=%s"
                 % (args.outBaseName,configDict['readFileName'],configDict['exonTabixFile']))

    pickreads.checkfile(configDict['readFileName'])
    pickreads.checkfile(configDict['exonTabixFile'])
    pickreads.checkfile(args.pgTabixFile)

    tabixFile = pysam.Tabixfile(configDict['exonTabixFile'], 'r')
    pgTabix = pysam.Tabixfile(args.pgTabixFile, 'r')

    # structure for keeping element classes seperate
    geneNameDict = {}

    f = open(configDict['readFileName'])
    for line in f:
        fields = line.rsplit("\t")
        chrom      = fields[0]
        readPos    = int(fields[1])
        readStrand = fields[2]
        mateElt    = fields[3]
        mateChrom  = fields[4]
        mateStrand = fields[5]
        matePos    = int(fields[6])
        eltStart   = int(fields[7])
        eltEnd     = int(fields[8])
        eltStrand  = fields[9]
        eltFullLen = int(fields[10])
        genomeName = fields[11]
        peakIndex  = fields[12]
        geneName   = fields[13].strip()

        if not geneNameDict.has_key(geneName):
            geneNameDict[geneName] = PeakBuilder(geneName);

        pair = pairinfo.Pair(chrom,readPos,readStrand,mateElt,mateChrom,
                             mateStrand,matePos,eltStart,eltEnd,eltStrand,
                             eltFullLen,genomeName,peakIndex,geneName)
        if geneNameDict[geneName].lastPeakNum != peakIndex:
            geneNameDict[geneName].peaks.append(geneNameDict[geneName].peak)
            geneNameDict[geneName].peak = pairinfo.Peak()
            geneNameDict[geneName].peak.addpair(pair)
            geneNameDict[geneName].lastPeakNum = peakIndex
        else:
            geneNameDict[geneName].peak.addpair(pair)

    for geneName in geneNameDict.keys():
        geneNameDict[geneName].peaks.append(geneNameDict[geneName].peak)

    logging.info("starting long output...")
    longOutput(geneNameDict,args.outBaseName,args.outDirName,pgTabix)

if __name__ == '__main__':
    # commandline args 
    parser = argparse.ArgumentParser(description='parse the output of pickreads.py')
    parser.add_argument('-c', '--config', dest='configFileName', default='config.txt',
                        help='config file left by pickreads.py')
    parser.add_argument('-o', '--outbasename', dest='outBaseName', required=True,
                        help='basename for output files')
    parser.add_argument('-d', '--outdirname', dest='outDirName', default='output',
                        help='output directory')
    parser.add_argument('-p', '--pseudogene', dest='pgTabixFile', required=True,
                        help='tabix file with pseudogene locations')
    args = parser.parse_args()

    main(args)
