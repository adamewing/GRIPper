#!/usr/bin/env python

"""
 Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)

 Released under the MIT license, see LICENSE.txt
"""


import logging
import os
import re
import argparse
import peakparser
import pickreads 

def cmpChrPosList(a,b):
    """
    +1 if a > b
    -1 if a < b
    0 if a == b
    """
    chrA = a[0]
    chrB = b[0]

    posA = int(a[1])
    posB = int(b[1])

    if (chrA == chrB):
        if (posA == posB): return 0
        if (posA < posB): return -1
        if (posA > posB): return 1

    if (chrA < chrB):
        return -1
    if (chrA > chrB):
        return 1

def sortChrPos(file):
    """Sorts files whose first two columns are chromosome and position, respectively."""
    f = open(file,'r')
    bedlist = []
    for bedline in f:
        bedlist.append(bedline.rstrip().split("\t"))
    return sorted(bedlist, cmp=cmpChrPosList)

def mergeChrPosFiles(infile1,infile2,outfile,maxDist):
    """Concatenates files, sorts, and outputs to outfile"""
    bedlist = []
    for file in (infile1,infile2):
        f = open(file,'r')
        findex = 0
        bedhead = ''
        for bedline in f:
            if findex == 0 and re.search(".bed$", outfile) and re.search("track", bedline): # ignore BED header
                bedhead = bedline
            else:
                bedlist.append(bedline.rstrip().split("\t"))
            findex = findex + 1
        f.close()

    outlist = []

    if re.search(".readpairs.txt$", outfile):
        outlist = reindexReads(sorted(bedlist, cmp=cmpChrPosList), maxDist)
    else:
        outlist = sorted(bedlist, cmp=cmpChrPosList)

    f = open(outfile, 'w')
    if re.search("track", bedhead):
        f.write("track\tname=%s\tvisibility=2\titemRgb=On\tdb=hg19\n" % os.path.basename(outfile))
    for outline in outlist: #sorted(bedlist, cmp=cmpChrPosList):
        outstring = "\t".join(outline)
        f.write(outstring + "\n")
    f.close()

class clusterBuild:
    """Stores placeholder data for building clusters"""
    def __init__(self):
        self.lastPos   = 0
        self.lastChr   = ''
        self.peakIndex = 0

def reindexReads(readlist, maxDist):
    """Used for generating new index numbers for peaks called by pickreads.py"""
    newlist = []

    cluster = clusterBuild();

    for readline in readlist:
        (readChrName,readPos,readStrand,mateElt,mateChrom,
         mateStrand,matePos,eltStart,eltEnd,eltStrand,
         eltFullLen,source,oldIndex,eltClass) = readline

        readPos = int(readPos)

        if (readChrName == cluster.lastChr):
            if (readPos - cluster.lastPos < 0): # sanity check
                raise IndexError('.bam file does not appear to be properly sorted')
            if (readPos - cluster.lastPos < maxDist): # keep the same peak number
                cluster.lastPos = readPos
            else: # create a new peak and add the pair to it if pair is > maxDist from last one
                cluster.peakIndex = cluster.peakIndex + 1
                cluster.lastPos = readPos
        else: # create a new peak and add the pair to it if pair is on a new chromosome
            cluster.peakIndex = cluster.peakIndex + 1
            cluster.lastChr = readChrName
            cluster.lastPos = readPos

        newlist.append((readChrName,str(readPos),readStrand,mateElt,mateChrom,
                        mateStrand,matePos,eltStart,eltEnd,eltStrand,eltFullLen,
                        source,str(cluster.peakIndex),eltClass))
    return newlist


def main(args):
    # debugging info
    logfile = args.outDirName + "/" + args.outBaseName + "/logs/%d" % os.getpid() + "." + args.outBaseName + ".mergepairs.log"
    logging.basicConfig(format='%(asctime)s %(message)s',filename=logfile,level=logging.DEBUG)

    logging.info("\ninDir1=%s\ninDir2=%s\noutBaseName=%s\nconfigFileName=%s"
                 % (args.inDir1,args.inDir2,args.outBaseName,args.configFileName))

    # create output directory
    pickreads.prepOutDir(args.outBaseName,args.outDirName,args.overwrite)

    # make sure input sources exist
    peakparser.checkOutDir(args.inDir1,args.outDirName)
    peakparser.checkOutDir(args.inDir2,args.outDirName)

    # make sure config files exist
    configPath1 = args.outDirName + "/" + args.inDir1 + "/" + args.configFileName
    configPath2 = args.outDirName + "/" + args.inDir2 + "/" + args.configFileName
    pickreads.checkfile(configPath1)
    pickreads.checkfile(configPath2)

    # read parameters for both inputs
    configDict1 = peakparser.readConfig(configPath1,args.inDir1,args.outDirName)
    configDict2 = peakparser.readConfig(configPath2,args.inDir2,args.outDirName)

    maxDist = int(configDict1['insertSize']) + 2*int(configDict1['readLength'])

    # merge readfiles
    outReadFileName = args.outDirName + "/" + args.outBaseName + "/" + args.outBaseName + ".readpairs.txt"
    readFileName1 = args.outDirName + "/" + args.inDir1 + "/" + args.inDir1 + ".readpairs.txt"
    readFileName2 = args.outDirName + "/" + args.inDir2 + "/" + args.inDir2 + ".readpairs.txt"
    logging.info("merging readfiles (%s, %s)" % (readFileName1,readFileName2))
#   print "merging readfiles (%s, %s)" % (readFileName1,readFileName2)
    mergeChrPosFiles(readFileName1,readFileName2,outReadFileName,maxDist)

    # merge bedfiles
    outBedFileName = args.outDirName +  "/" + args.outBaseName + "/" + args.outBaseName + ".reads.bed"
    bedFileName1 = args.outDirName + "/" + args.inDir1 + "/" + args.inDir1 + ".reads.bed"
    bedFileName2 = args.outDirName + "/" + args.inDir2 + "/" + args.inDir2 + ".reads.bed"
    logging.info("merging bedfiles (%s,%s)" % (bedFileName1,bedFileName2))
#   print "merging bedfiles (%s,%s)" % (bedFileName1,bedFileName2)
    mergeChrPosFiles(bedFileName1,bedFileName2,outBedFileName,maxDist)

    # write new config file
    configPath = args.outDirName + "/" + args.outBaseName + "/" + args.configFileName
    configDict = configDict1
    configDict['bamFileName1'] = configDict1['bamFileName']
    configDict['bamFileName2'] = configDict2['bamFileName']
    configDict['merged'] = 'True'
    configDict['outBaseName'] = args.outBaseName
    configDict['outDirName'] = args.outDirName
    configDict['readFileName'] = outReadFileName

    del configDict['bamFileName']

    f = open(configPath, 'w')
    for k,v in configDict.iteritems():
        f.write(k + "=" + v + "\n")
    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='merge files in the specified output directories, maintains sorting')
    parser.add_argument('-1', '--1', dest='inDir1', required=True,
                        help='output directory 1 (input 1)')
    parser.add_argument('-2', '--2', dest='inDir2', required=True,
                        help='output directory 2 (input 2)')
    parser.add_argument('-o', '--outbasename', dest='outBaseName', required=True,
                        help='basename for merged output')
    parser.add_argument('-d', '--outdirname', dest='outDirName', default='output',
                        help='output directory')
    parser.add_argument('-c', '--config', dest='configFileName', default='config.txt',
                        help='filname for config file left by discordant.py')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()

    main(args)
