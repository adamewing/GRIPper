#!/usr/bin/env python

"""
 Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)

 Released under the MIT license, see LICENSE.txt
"""

import pysam
import argparse
import sys
import logging
import os
import re

def checkfile(fname):
    """check whether a file exists"""
    try:
        open(fname)
    except IOError as e:
        print "can't find file: " + fname
        sys.exit()

def prepOutDir(outBaseName,outDirName,overwrite):
    if not os.path.exists(outDirName):
        os.mkdir(outDirName)
    if not os.path.exists(outDirName + "/" + outBaseName):
        os.mkdir(outDirName + "/" + outBaseName)
    else:
        if not overwrite:
            sys.exit(outDirName + " exists, can be overridden with --overwrite")
    if not os.path.exists(outDirName + "/" + outBaseName + "/logs"):
        os.mkdir(outDirName + "/" + outBaseName + "/logs")

def writeConfig(args):
    """Write configuration information to a config.txt"""
    f = open(args.outDirName + "/" + args.outBaseName + "/config.txt", 'w')
    for k,v in vars(args).iteritems():
         f.write(k + "=" + str(v) + "\n")
    f.close()

def fixChrName(name):
    """fix for the most common variation on chromosome names (leaving out the 'chr')"""
    if name[0] != 'c':
        return "chr" + name
    else:
        return name

class clusterBuild:
    """Stores placeholder data for building clusters"""
    def __init__(self):
        self.lastPos   = 0
        self.lastChr   = ''
        self.peakIndex = 0

def main(args):
    prepOutDir(args.outBaseName,args.outDirName,args.overwrite)
    writeConfig(args) # save arguments in output directory

    # debugging info
    logfile = args.outDirName + "/" + args.outBaseName + "/logs/%d" % os.getpid() + "." + args.outBaseName +".discordant.log"
    logging.basicConfig(format='%(asctime)s %(message)s',filename=logfile,level=logging.DEBUG)

    # log parameters
    logging.info("\noutBaseName=%s\nbamFileName=%s\nexonTabixFile=%s\ngeneTabixFile=%s\nreadLength=%s\ninsertSize=%s" 
             % (args.outBaseName,args.bamFileName,args.exonTabixFile,args.geneTabixFile,args.readLength,args.insertSize))

    checkfile(args.bamFileName)
    checkfile(args.exonTabixFile)
    checkfile(args.geneTabixFile)

    bamFile = pysam.Samfile(args.bamFileName, 'rb') # rb = read, binary
    tabixFile = pysam.Tabixfile(args.exonTabixFile, 'r')
    geneTabixFile = pysam.Tabixfile(args.geneTabixFile, 'r')
    tabixContigs = {}
    for contig in tabixFile.contigs:
        tabixContigs[contig] = 1

    # parameters for reads
    readLength = int(args.readLength)
    insertSize = int(args.insertSize)
    maxDist = readLength + (2*insertSize)

    bamIndex = 0

    # output file
    readsFileName = args.outBaseName + ".readpairs.txt"
    bedFileName   = args.outBaseName + ".reads.bed"
    readsFile = open(args.outDirName + "/" + args.outBaseName + "/" + readsFileName, 'w')
    bedFile   = open(args.outDirName + "/" + args.outBaseName + "/" + bedFileName, 'w')

    bedHeader = "track\tname=" + args.outBaseName + "_Reads\tvisibility=2\titemRgb=On\tdb=" + args.refGenome + "\n"
    bedFile.write(bedHeader)

    cluster = clusterBuild()

    for read in bamFile.fetch():
    # avoid putting anything else before this 'if statement' as it will be evaluated for _every_ read
    # in the .bam file (and there are a lot of reads...)

        # we'll get the mates of reads that are mapped to insertion site coords since pysam can't return the sequence of the mate
        if (not read.is_proper_pair and not read.is_duplicate and not read.is_unmapped 
            and read.mapq >= 20 and read.tid > -1 and read.mrnm > -1 
            and tabixContigs.has_key(fixChrName(bamFile.getrname(read.tid)))
            and tabixContigs.has_key(fixChrName(bamFile.getrname(read.mrnm)))):

        # eventually, we'll want to keep track of two clusters:
        # 1. exon-exon junctions that are discordant relative to the reference
        # 2. anchors with discordant/unmapped mates

            readChrName = fixChrName(bamFile.getrname(read.tid))
            mateChrName = fixChrName(bamFile.getrname(read.mrnm))
            seqLen = len(read.seq)

            # fetch() using the pysam.asTuple parser returns tabix results as a tuple
            tabixTupleParse = tabixFile.fetch(reference=readChrName, 
                                              start=read.pos, 
                                              end=read.pos+seqLen, 
                                              parser=pysam.asTuple())

            tabixMateTupleParse = tabixFile.fetch(reference=mateChrName, 
                                                  start=read.mpos, 
                                                  end=read.mpos+seqLen, 
                                                  parser=pysam.asTuple())

            readExon = None 
            mateExon = None

            if tabixTupleParse:
                for tabixTuple in tabixTupleParse:
                    readExon = tabixTuple[3]
                for tabixMateTuple in tabixMateTupleParse:
                    mateExon = tabixMateTuple[3]

            readStrand = '+'
            if read.is_reverse:
                readStrand = '-'

            mateStrand = '+'
            if read.mate_is_reverse:
                mateStrand = '-'

            if (readChrName == cluster.lastChr):
                if (read.pos - cluster.lastPos < 0): # sanity check
                    raise IndexError('.bam file does not appear to be properly sorted')
                if (read.pos - cluster.lastPos < maxDist): # keep the same peak number
                    cluster.lastPos = read.pos
                else: # create a new peak and add the pair to it if pair is > maxDist from last one
                    cluster.peakIndex = cluster.peakIndex + 1
                    cluster.lastPos = read.pos
            else: # create a new peak and add the pair to it if pair is on a new chromosome
                cluster.peakIndex = cluster.peakIndex + 1
                cluster.lastChr = readChrName
                cluster.lastPos = read.pos

            if readExon == None and mateExon != None:
                # get information about gene
                geneStrand = ''
                geneMaxLen = 0

                geneTabixParse = geneTabixFile.fetch(reference=mateChrName,
                                                     start=read.mpos,
                                                     end=read.mpos+seqLen,
                                                     parser=pysam.asTuple())

                if geneTabixParse:
                    for geneInfo in geneTabixParse:
                        glen = int(geneInfo[2])-int(geneInfo[1])
                        if glen > geneMaxLen:
                            geneMaxLen = glen
                            geneStrand = geneInfo[4]
                            geneName   = geneInfo[3]

                out = "\t".join((readChrName,str(read.pos),readStrand,str(mateExon),
                                 mateChrName,mateStrand,str(read.mpos),tabixMateTuple[1],tabixMateTuple[2],
                                 geneStrand,str(geneMaxLen),args.outBaseName,str(cluster.peakIndex),geneName))
                readsFile.write(out + "\n")

        bamIndex = bamIndex + 1
    readsFile.close()
    bedFile.close()
    logging.info("%s finished, (%d reads)" % (args.bamFileName,bamIndex))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='parse .bam file and output discordant reads where one end is in an annotation specified in a tabix index')
    parser.add_argument('-b', '--bam', dest='bamFileName', required=True, 
                        help='name of .bam file')
    parser.add_argument('-t', '--tabix', dest='exonTabixFile', required=True,
                        help='name of exon tabix index (.gz)') 
    parser.add_argument('-e', '--genetabix', dest='geneTabixFile', required=True,
                        help='name of gene tabix index')
    parser.add_argument('-o', '--outbasename', dest='outBaseName', default='defaultName',
                        help='basename for output files')
    parser.add_argument('-d', '--outdirname', dest='outDirName', default='output',
                        help='output directory')
    parser.add_argument('-r', '--readlen', dest='readLength', default='100',
                        help='read length in basepairs (default  100 bp)')
    parser.add_argument('-i', '--insertsize', dest='insertSize', default='300',
                        help='expected insert size in bp (default 300 bp)')
    parser.add_argument('-g', '--refgenome', dest='refGenome', default='hg19',
                        help='reference genome assembly (defaults to hg19)')
    parser.add_argument('--overwrite', action="store_true")
    args = parser.parse_args()

    main(args)
