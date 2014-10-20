#!/usr/bin/env python

"""
 Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)

 Released under the MIT license, see LICENSE.txt
"""

import sys, re, time, pp, argparse, subprocess, pysam, ConfigParser
import lib.configtest as configtest
from os import path as path

# other imports are done via pp

class PairedSample:
    def __init__(self,name,config,outdir,assembly,tcga=False,justcall=False):
        self.name = name
        self.tcga = tcga      # if true, changes validation, naming
        self.justcall = justcall
        self.config = config  # ConfigParser
        self.outdir = outdir
        self.cancerBam = None
        self.normalBam = None
        self.cancerIdx = None
        self.normalIdx = None 
        self.assembly  = assembly 

    def __str__(self):
        return "\t".join((str(self.cancerBam),str(self.normalBam),
                          str(self.cancerIdx),str(self.normalIdx)))

    def validate(self):
        if self.tcga:
            if (self.cancerBam and path.exists(self.cancerBam) and
                self.normalBam and path.exists(self.normalBam) and
                self.cancerIdx and path.exists(self.cancerIdx) and
                self.normalIdx and path.exists(self.normalIdx)):
                return True
            else:
                return False
        else:
            if (self.normalBam and path.exists(self.normalBam) and
                self.normalIdx and path.exists(self.normalIdx)):
                return True
            else:
                return False

    def addFile(self,sampleType,extension,filePath):
        if sampleType == 'NORMAL' and extension == 'bam':
            self.normalBam = filePath
            self.normalIdx = filePath + ".bai"
        if sampleType == 'CANCER' and extension == 'bam':
            self.cancerBam = filePath
            self.cancerIdx = filePath + ".bai"


    def runDiscordant(self):
        """Runs pickreads.py, peakparser.py, summarize,py, pinpoint.py"""
        exonTabix = self.config.get('discord', 'tabixDir') + "/" + self.assembly + ".exons.bed.gz"
        geneTabix = self.config.get('discord', 'tabixDir') + "/" + self.assembly + ".geneinfo.txt.gz"
        usechr = False
        if self.config.get('discord', 'usechr') == "True":
            usechr = True
        args = argparse.Namespace(bamFileName    = self.normalBam,  # pickreads.py
                                  outBaseName    = self.name,       # pickreads.py, peakparser.py, summarize.py, pinpoint.py
                                  exonTabixFile  = exonTabix,       # pickreads.py
                                  geneTabixFile  = geneTabix,       # pickreads.py
                                  refGenome      = self.assembly,   # pickreads.py
                                  overwrite      = True,            # pickreads.py
                                  printout       = False,           # summarize.py
                                  inDir1         = None,            # mergepairs.py
                                  inDir2         = None,            # mergepairs.py
                                  outDirName     = self.outdir,     # everything
                                  usechr         = usechr,          # pinpoint.py
                                  readLength     = self.config.get('discord','readLength'),
                                  insertSize     = self.config.get('discord','insertSize'),
                                  minMapQ        = self.config.get('discord','minMapQ'),
                                  configFileName = self.config.get('discord','configFileName'),
                                  minPeakSize    = self.config.get('discord','minPeakSize'),
                                  maxReadLen     = self.config.get('discord','maxReadLen'),
                                  zeroChar       = self.config.get('discord','zeroChar'),
                                  minClipQual    = self.config.get('discord','minClipQual'),
                                  mrnaFastaFile  = self.config.get('discord','mrnaFastaFile'),
                                  annotDir       = self.config.get('discord','annotDir'),
                                  refGenomeFile  = self.config.get('discord',self.assembly),
                                  pgTabixFile    = self.config.get('discord','pgTabixFile'))

        if self.tcga:
            normalName = self.name + "-NORMAL"
            cancerName = self.name + "-CANCER"
            mergeName  = self.name + "-MERGE"
            args.inDir1 = normalName
            args.inDir2 = cancerName
            args.outBaseName = normalName
    
            # normal .bam
            lib.pickreads.main(args)

            # cancer .bam
            args.bamFileName = self.cancerBam
            args.outBaseName = cancerName
            lib.pickreads.main(args)

            # merged results
            args.outBaseName = mergeName
            lib.mergepairs.main(args)

            lib.peakparser.main(args)
            lib.summarize.main(args)
            lib.pinpoint.main(args)

        else:
            if not self.justcall:
                lib.pickreads.main(args)
            lib.peakparser.main(args)
            lib.summarize.main(args)
            lib.pinpoint.main(args)

def main(args):
    configtest.check(args.configFile)
    
    pairedSamples = {}

    for line in open(args.sampleFile, 'r'):
        if not re.search("^#", line):
            (filePath,label,assembly) = line.strip().split()
            base = path.basename(filePath)

            sampleType = 'NORMAL'
            sampleName = label
            fileparts = base.split('.')
            extension = fileparts[-1]

            cancerlabellist = ['tumor', 'tumour', 'cancer']
            for cancerlabel in cancerlabellist:
                if re.search(cancerlabel, sampleName.lower()):
                    sampleType = 'CANCER'

            if args.tcga: # see https://wiki.nci.nih.gov/display/TCGA/Working+with+TCGA+Data
                baseparts = base.split('-')
                participant = baseparts[2]
                samplenum = baseparts[3]
                if re.search('01', samplenum):
                    sampleType = 'CANCER'
                sampleName = '-'.join((label,participant))

            config = ConfigParser.ConfigParser()
            config.read(args.configFile)

            if sampleName not in pairedSamples:
                pairedSamples[sampleName] = PairedSample(sampleName,config,args.outDirName,assembly,tcga=args.tcga,justcall=args.justcall)

            pairedSamples[sampleName].addFile(sampleType,extension,filePath)

    # parallel python stuff
    ncpus = int(args.numCPUs) 
    jobServer = pp.Server(ncpus,ppservers=())

    sampleJobs = {}
    for sampleName in pairedSamples.keys():
        if pairedSamples[sampleName].validate():
            print sampleName + " validated"
            job = jobServer.submit(pairedSamples[sampleName].runDiscordant,(),(),
                  ("lib.pickreads","lib.mergepairs","lib.peakparser","lib.pinpoint","lib.summarize","argparse"))
            sampleJobs[sampleName] = job
            print "started job: " + sampleName
            jobServer.print_stats()
        else:
            print sampleName + " did not validate, make sure the .bam exists and is indexed (has a .bam and a .bam.bai file)"

    for sampleName in sampleJobs:
        print "=== stdout of job: " + sampleName + " ==="
        print sampleJobs[sampleName]()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Pipeline that aims to identify TE insertions from .bam alignments")
    parser.add_argument('-s', '--sampleFile', dest="sampleFile", required=True,
                    help="List of filenames and sample names. Should include both .bam files (.bam) and their indexes (.bai)")
    parser.add_argument('-c', '--config', dest='configFile', required=True)
    parser.add_argument('-p', '--processors', dest='numCPUs', default='1', help="Number of CPUs to use (1 job per CPU)")
    parser.add_argument('-o', '--outdir', dest='outDirName', required=True, help="output base directory")
    parser.add_argument('--tcga', action="store_true", help="Use TCGA filename format to determine cancer/normal pair")
    parser.add_argument('--justcall', action="store_true", help="don't pick discordant reads, just call on pre-existing sample.readpairs.txt")
    args = parser.parse_args()
    main(args)
