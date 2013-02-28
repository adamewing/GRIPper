#!/usr/bin/env python

"""
 Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)

 Released under the MIT license, see LICENSE.txt
"""

import re
import os
import subprocess
import tempfile
import pysam
import sys
import gzip
import argparse
import peakparser
from string import maketrans

def checkfile(fname):
    try:
        open(fname)
    except IOError as e:
        print "can't find file: " + fname
        sys.exit()

# parse TCGA filename
def getTypeFromTCGA(fname):
    try:
        fields = os.path.basename(fname).split('-')
        sample = fields[3]
        if sample[0] == '0':
            return 'CANCER'
        if sample[0] == '1':
            return 'NORMAL'
        return None
    except:
        #print "invalid TCGA filename: " + fname
        return None

# some people choose to use reference genomes with chromosome names that don't begin in 'chr'
# since we do use 'chr', we need to detect those cases
def chromNameUsesPrefix(bam):
    for ref in bam.references:
        if re.search("^chr",ref):
            return True
    return False

# fix for the most common variation on chromosome names (leaving out the 'chr')
def fixChrName(name):
    if name[0] != 'c':
        return "chr" + name
    else:
        return name

# average base quality over an interval
def avgQual(qstring,start,end,zeroChar):
    if start < 0:
        start = 0
    if end > len(qstring):
        end = len(qstring)

    offset = ord(zeroChar)

    return sum(map(lambda x: ord(x)-offset, qstring[start:end]))/(end-start)

# capitalize 1-indexed seq between start and end
def capSeq(seq,start,end):
    start -= 1
    seq = seq.lower()
    chars = []
    for c in seq:
        chars.append(c)
    for i in range(start,end):
        chars[i] = chars[i].upper()
    return ''.join(chars)

# runs bwa stdsw to align two sequences
def bwastdsw(queryName,query,refName,ref,queryIsFile=False,refIsFile=False):
    qfileName = ''
    rfileName = ''

    if queryIsFile:
        qfileName = query
    else:
        qfile = tempfile.NamedTemporaryFile(delete=False)
        qfile.write(">%s\n%s\n" % (queryName,query))
        qfile.close()
        qfileName = qfile.name

    if refIsFile:
        rfileName = ref
    else:
        rfile = tempfile.NamedTemporaryFile(delete=False)
        rfile.write(">%s\n%s\n" % (refName,ref))
        rfile.close()
        rfileName = rfile.name

    # can use the -f option to only consider forward strand since bwa prints aligned seqs
    # in the forward direction
    args = ['bwa', 'stdsw', rfileName, qfileName]
    p = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)

    fnum = 0
    alignments = []
    for pline in p.stdout.readlines():
        if re.search("^>", pline):
            fnum = 0
            col = pline.strip().split("\t")
            alignments.append(BwastdswAlignResult(col))
            alignments[-1].queryLen = len(query) # FIXME make this work with filename input
            alignments[-1].targetLen = len(ref)
           # print "%d\t" % fnum + "\t".join(col)
            fnum += 1
        elif fnum == 1:
            alignments[-1].targetSeq = pline.strip()
           # print "%d\t" % fnum + pline.strip()
            fnum += 1
        elif fnum == 2:
            alignments[-1].matchString = pline.strip()
           # print "%d\t" % fnum + pline.strip()
            fnum += 1
        elif fnum == 3:
            alignments[-1].querySeq = pline.strip()
           # print "%d\t" % fnum + pline.strip()
            fnum += 1
    p.stdout.close()
    p.kill()
    p.wait()

    if not queryIsFile:
        os.unlink(qfileName)
    if not refIsFile:
        os.unlink(rfileName)

    if len(alignments) == 0:
        return None
    if len(alignments) == 1:
        return alignments[0]

    # more than 1 alignment, pick the 'best'
    bestaln = alignments[0] 
    for align in alignments:
        if align.alnLength() > bestaln.alnLength() and align.pctID() > 90:
            bestaln = align

    return bestaln

def fetchRegion(bamFile,refGenome,maxReadLen,chr,start,end,gene,zeroChar,minClipQual,usechr=False):
    maxReadLen = int(maxReadLen)
    start = int(start)
    end = int(end)

    if not usechr or not chromNameUsesPrefix(bamFile):
        chr = chr.lstrip('chr')

    regionRefSeq = refGenome.fetch(reference=chr, start=start-maxReadLen, end=end+maxReadLen)

    cluster = ClippedCluster(chr,start,end,gene,maxReadLen)

    for read in bamFile.fetch(chr, start-maxReadLen, end+maxReadLen):
        if not read.is_unmapped and not read.is_duplicate:
            cliplen = read.rlen - read.qlen # minimum soft clipping: if rlen < qlen --> bases were soft-clipped
            if cliplen > 10:
                refseq = ''
                leftclip = 0
                rightclip = 0
                if (read.qstart == cliplen): # left-clipped
                    leftclip = cliplen

                if (read.qstart == 0): # right-clipped
                    rightclip = cliplen

                if (read.qstart > 0 and read.qstart < cliplen):
                    leftclip = cliplen - read.qstart
                    rightclip = cliplen - leftclip

                breakside = ''
                breakloc  = 0
                clipqual  = 0

                if (leftclip > rightclip-10): # 10 is arbitrary
                    breakside = 'L'
                    breakloc = read.pos
                    clipqual = avgQual(read.qual,0,leftclip,zeroChar)
                elif (rightclip > leftclip-10):
                    breakside = 'R'
                    breakloc = read.pos + (read.rlen - rightclip)
                    clipqual = avgQual(read.qual,rightclip,read.rlen,zeroChar)
                else:
                    breakside = 'A' # ambiguous

                if (clipqual >= int(minClipQual)):
                    align = bwastdsw('query',read.seq,'target',regionRefSeq)
                    #print align
                    if align:
                        breakLeft  = start - maxReadLen + align.targetStart
                        breakRight = start - maxReadLen + align.targetEnd

                        if (breakLeft >= start-10 and breakLeft <= end+10) or (breakRight >= start-10 and breakRight <= end+10):
                            cluster.aligns.append(align)
                            cluster.reads.append(read)
    cluster.assignBreaks()
    return cluster;

def mergeClusters(cl1,cl2,txList):

    if not cl1.hasReads():
        cl2.mapTx(txList)
        return cl2

    if not cl2.hasReads():
        cl1.mapTx(txList)
        return cl1

    if (cl1.chr == cl2.chr and cl1.start == cl2.start and cl1.end == cl2.end and cl1.gene == cl2.gene and cl1.maxrdln == cl2.maxrdln):
        new = ClippedCluster(cl1.chr, cl1.start, cl1.end, cl1.gene, cl1.maxrdln)

        new.aligns = cl1.aligns
        for align in cl2.aligns:
            new.aligns.append(align)

        new.reads = cl1.reads
        for read in cl2.reads:
            new.reads.append(read)

        new.assignBreaks()
        new.mapTx(txList)
        new.type = cl1.type + "," + cl2.type
        return new
    else:
        raise ValueError('cannot merge clusters that have different locations')

class ClippedCluster:
    def __init__(self,chr,start,end,gene,maxReadLen):
        self.chr     = chr
        self.start   = int(start)
        self.end     = int(end)
        self.gene    = gene
        self.maxrdln = maxReadLen
        self.aligns  = [] # BwastdswAlignResult objects
        self.reads   = [] # pysam AlignedRead objects
        # assigned by functions:
        self.assign  = False # has assignBreaks() been run?
        self.lgood   = False # is left break good? (defined in bestBreakLeft())
        self.rgood   = False # is right break good? (defined in bestBreakRight())
        self.lbest   = 0 # best guess for left break
        self.rbest   = 0 # best guess for right break
        self.lbreaks = [] # candidate left break positions
        self.rbreaks = [] # candidate right break positions
        self.type = ''

        self.teqlen  = [] # TE query lengths
        self.tetype  = [] # TE families (list of those detected)
        self.testart = [] # list of TE starts (positions in TEs)
        self.teend   = [] # list of TE ends (positions in TEs)
        self.testr   = [] # list of TE orientations
        self.minTEQueryLen = 15 # minimum length of TE seq for alignment to be valid


    def hasReads(self):
        if len(self.reads) > 0:
            return True
        return False

    def mapTx(self,txs):
        for i in range(len(self.reads)):
            teAlign = partialMapTx(txs,
                                   self.gene,
                                   self.reads[i].seq,
                                   self.aligns[i].queryStart,
                                   self.aligns[i].queryEnd)

            if (teAlign == None):
                self.tetype.append('None')
                self.testart.append(0)
                self.teend.append(0)
                self.testr.append('.')
                self.teqlen.append(0)
            elif (teAlign.queryLen >= self.minTEQueryLen):
                self.tetype.append(teAlign.targetName)
                self.testart.append(teAlign.targetStart)
                self.teend.append(teAlign.targetEnd)
                self.testr.append(teAlign.queryStr)
                self.teqlen.append(teAlign.queryLen)
            else:
                self.tetype.append('None')
                self.testart.append(0)
                self.teend.append(0)
                self.testr.append('.')
                self.teqlen.append(teAlign.queryLen)

    def assignBreaks(self):
        for i in range(len(self.reads)):
            leftmargin  = self.aligns[i].queryStart - 1
            rightmargin = self.reads[i].rlen - self.aligns[i].queryEnd
            if leftmargin > rightmargin and rightmargin <= 10: 
                self.lbreaks.append(self.start - self.maxrdln + self.aligns[i].targetStart) 
            if leftmargin < rightmargin and leftmargin <= 10:
                self.rbreaks.append(self.start - self.maxrdln + self.aligns[i].targetEnd) 
        self.assign = True

        if len(self.lbreaks) > 0:
            self.bestBreakLeft()
        if len(self.rbreaks) > 0:
            self.bestBreakRight()
        self.retryBreaks()

    def bestBreakLeft(self):
        if not self.assign:
            self.assignBreaks()
        # count unique positions
        posCountDict = {}
        for pos in self.lbreaks:
            if posCountDict.has_key(pos):
                posCountDict[pos] += 1
            else:
                posCountDict[pos] = 1
        # get modal position
        maxCount = 0
        maxCountPos = 0
        for pos,count in posCountDict.iteritems():
            if count > maxCount:
                maxCount = count
                maxCountPos = pos
        # if over half the reads have this break, it's good
        if float(maxCount)/float(len(self.lbreaks)) > 0.5:
            self.lgood = True
        self.lbest = int(maxCountPos)

    def bestBreakRight(self):
        if not self.assign:
            self.assignBreaks()
        # count unique positions
        posCountDict = {}
        for pos in self.rbreaks:
            if posCountDict.has_key(pos):
                posCountDict[pos] += 1
            else:
                posCountDict[pos] = 1
        # get modal position
        maxCount = 0
        maxCountPos = 0
        for pos,count in posCountDict.iteritems():
            if count > maxCount:
                maxCount = count
                maxCountPos = pos
        # if over half the reads have this break, it's good
        if float(maxCount)/float(len(self.rbreaks)) > 0.5:
            self.rgood = True
        self.rbest = int(maxCountPos)

    def retryBreaks(self):
        """look for breaks that aren't majority but do indicate a TSD"""

        if self.lgood and not self.rgood:
            # look for TSD in rbreaks
            tsdopts = {} # TSD candidates, stores count
            for rpos in self.rbreaks:
                if (rpos - self.lgood) >= 2 and (rpos - self.lgood) <= 50:
                    if rpos in tsdopts:
                        tsdopts[rpos] += 1
                    else:
                        tsdopts[rpos] = 1
            maxCount = 0
            maxCountPos = 0
            numTiedForBest = 0
            for pos,count in tsdopts.iteritems():
                if count > maxCount:
                    maxCount = count
                    maxCountPos = pos
                if count == maxCount:
                    numTiedForBest += 1
            if maxCount > 0 and numTiedForBest == 0:
                self.rgood = True
                sys.stderr.write("better break found for " + self.chr + ":" + str(self.start) + "-" + str(self.end) + "\n")
                self.rbest = int(maxCountPos)

        if self.rgood and not self.lgood:
            # look for TSD in lbreaks
            tsdopts = {} # TSD candidates, stores count
            for lpos in self.lbreaks:
                if (self.rgood - lpos) >= 2 and (self.rgood - lpos) <= 50:
                    if lpos in tsdopts:
                        tsdopts[lpos] += 1
                    else:
                        tsdopts[lpos] = 1
            maxCount = 0
            maxCountPos = 0
            numTiedForBest = 0
            for pos,count in tsdopts.iteritems():
                if count > maxCount:
                    maxCount = count
                    maxCountPos = pos
                if count == maxCount:
                    numTiedForBest += 1
            if maxCount > 0 and numTiedForBest == 0:
                self.lgood = True
                sys.stderr.write("better break found for " + self.chr + ":" + str(self.start) + "-" + str(self.end) + "\n")
                self.lbest = int(maxCountPos)

    def majorityTE(self):
        """returns the most frequently identified TE"""
        tecount = {}
        i = 0
        for te in self.tetype:
            # only count reads with long enough TE seq for alignment
            if te in tecount and self.teqlen[i] >= self.minTEQueryLen:
                tecount[te] += 1
            elif self.teqlen[i] >= self.minTEQueryLen:
                tecount[te] = 1
            i += 1

        if len(self.tetype) > 0:
            majTE = self.tetype[0]
            maxcount = 0
            for (te,count) in tecount.iteritems():
                if count > maxcount:
                    majTE = te
                    maxcount = count
            return majTE
        else:
            return 'None'

    def outstring(self):
        output = "\t".join((str(self.lgood), str(self.rgood), str(self.lbest), str(self.rbest), 
                            str(len(self.lbreaks)), str(len(self.rbreaks)), str(len(self.reads)), 
                            self.type, self.majorityTE()))
        return output

    def infodump(self):
        output = ''
        for i in range(len(self.reads)):
            rbreak = self.start - self.maxrdln + self.aligns[i].targetEnd
            lbreak = self.start - self.maxrdln + self.aligns[i].targetStart
            outseq = capSeq(self.reads[i].seq, self.aligns[i].queryStart, self.aligns[i].queryEnd)
            output += ("%s tr=%d tl=%d lr=%d ll=%d te=%s,%s,%s,%s"
                    % (outseq, rbreak, lbreak, self.aligns[i].queryStart, self.aligns[i].queryEnd,
                       self.tetype[i],str(self.testart[i]),str(self.teend[i]),self.testr[i]) + "\n")
        output += ("leftgood=%s rightgood=%s leftbreak=%d rightbreak=%d type=%s"
                % (self.lgood, self.rgood, self.lbest, self.rbest, self.type) + "\n")
        return output

class BwastdswAlignResult:
    def __init__(self,col):
        self.queryName   = col[3].lstrip('>')
        self.targetName  = col[0].lstrip('>')
        self.queryStart  = int(col[1])
        self.queryEnd    = int(col[2])
        self.queryStr    = col[4]
        self.targetStart = int(col[5])
        self.targetEnd   = int(col[6])
        self.cigarString = col[8]
        self.queryLen    = 0
        self.targetLen   = 0
        self.targetSeq   = ''
        self.querySeq    = ''
        self.matchString = ''

    def pctID(self):
        return float(sum(map(lambda x: int(x=='|'), self.matchString)))/float(len(self.querySeq))*100
    def alnLength(self):
        return len(self.targetSeq)
    def __str__(self):
        output = ("queryName=%s targetName=%s queryStart=%d queryEnd=%d queryStr=%s targetStart=%d targetEnd=%d pctID=%f queryLen=%d targetLen=%d"
                 % (self.queryName, self.targetName, self.queryStart, self.queryEnd, self.queryStr, self.targetStart, self.targetEnd, self.pctID(),
                 self.queryLen, self.targetLen))
        return output

class TranscriptSeq:
    def __init__(self,genename):
        self.gene = genename
        self.txs  = [] # transcript numbers
        self.seqs = []
        self.junc = []
    def __str__(self):
        output = self.gene + "\n"
        for i in range(len(self.txs)):
            output += str(self.txs[i]) + "\t" + self.seqs[i] + "\n"
        return output

def fastahash(infile):
    f = gzip.open(infile,'r')
    tx = {}
    gene = ''
    seq  = ''
    for line in f:
        if re.search("^>",line):
            if gene:
                tx[gene].seqs.append(seq) # finish last gene
            fields = line.strip().strip(">").split('.')
            juncstr = fields[-1]
            txnum = fields[-2]
            gene = fields[-3]
            if tx.has_key(gene):
                tx[gene].txs.append(int(txnum))
                tx[gene].junc.append(juncstr)
            else:
                tx[gene] = TranscriptSeq(gene)
                tx[gene].txs.append(int(txnum))
                tx[gene].junc.append(juncstr)
            seq = ''
        else:
            seq += line.strip()
    f.close()
    tx[gene].seqs.append(seq)
    return tx


def partialMapTx(txs,queryName,querySeq,excludeStart,excludeEnd):
    """
    Maps querySeq to a list of FASTA files (with names), returns best result (as BwastdswAlignResult)
    Only maps region outside of (excludeStart,excludeEnd)
    """
    excludeStart = int(excludeStart)
    excludeEnd = int(excludeEnd)
    results = []

    if queryName in txs:
        for txseq in txs[queryName].seqs:
            leftQuery  = querySeq[0:excludeStart]
            rightQuery = querySeq[excludeEnd:len(querySeq)-1]
            leftAlign  = bwastdsw(queryName,leftQuery,queryName,txseq)
            rightAlign = bwastdsw(queryName,rightQuery,queryName,txseq)
            if leftAlign != None and leftAlign.pctID > 90 and len(leftAlign.querySeq) >= 15:
                results.append(leftAlign)
            if rightAlign != None and rightAlign.pctID > 90 and len(rightAlign.querySeq) >= 15:
                results.append(rightAlign)

    bestPctID = 0
    bestResult = None
    for result in results:
        if result != None and result.pctID() > bestPctID:
            bestResult = result
            bestPctID = result.pctID()

    return bestResult

def main(args):
    peakparser.checkOutDir(args.outBaseName,args.outDirName)

    configPath = args.outDirName + "/" + args.outBaseName + "/" + args.configFileName
    checkfile(configPath)
    configDict = peakparser.readConfig(configPath,args.outBaseName,args.outDirName)

    # load hash of TranscriptSeq objects (genename --> multiple transcripts)
    # will need to pass to partialMapTx() later
    sys.stderr.write("loading " + args.mrnaFastaFile + "...\n")
    checkfile(args.mrnaFastaFile)
    txs = fastahash(args.mrnaFastaFile)

    cancerBamFile = ''
    normalBamFile = ''
    refGenomeFile = args.refGenomeFile 

    cancerCallsFile = args.outDirName + "/" + args.outBaseName + "/canceronly.tab.txt"
    normalCallsFile = args.outDirName + "/" + args.outBaseName + "/normalonly.tab.txt"
    germCallsFile   = args.outDirName + "/" + args.outBaseName + "/germline.tab.txt"
    otherCallsFile  = args.outDirName + "/" + args.outBaseName + "/uncategorized.tab.txt"

    # fix if unmerged
    if not configDict.has_key('bamFileName1'):
        configDict['bamFileName1'] = configDict['bamFileName']
        configDict['bamFileName2'] = configDict['bamFileName']

    bamType1 = getTypeFromTCGA(configDict['bamFileName1']) 
    bamType2 = getTypeFromTCGA(configDict['bamFileName2'])

    print ("bamfile1=%s bamFile2=%s bamType1=%s bamType2=%s"
        % (configDict['bamFileName1'], configDict['bamFileName2'], bamType1, bamType2))

    if bamType1 != bamType2 and bamType1 != None and bamType2 != None:
        if bamType1 == 'CANCER':
            if bamType2 != 'NORMAL':
                raise NameError('bam1 is cancer but bam2 is not normal')
            cancerBamFile = configDict['bamFileName1']
            normalBamFile = configDict['bamFileName2']
        if bamType2 == 'CANCER':
            if bamType1 != 'NORMAL':
                raise NameError('bam2 is cancer but bam1 is not normal')
            cancerBamFile = configDict['bamFileName2']
            normalBamFile = configDict['bamFileName1']
    else:
        print 'cannot determine bamfile cancer/normal from filenames in config.txt, defaulting to normal.'
        normalBamFile = configDict['bamFileName1']
        cancerBamFile = configDict['bamFileName2']

    checkfile(cancerBamFile)
    checkfile(normalBamFile)
    checkfile(normalCallsFile)
    checkfile(cancerCallsFile)
    checkfile(germCallsFile)
    checkfile(otherCallsFile)
    checkfile(refGenomeFile)

    cancerBam   = pysam.Samfile(cancerBamFile, 'rb') # rb = read, binary
    normalBam   = pysam.Samfile(normalBamFile, 'rb') # rb = read, binary
    cancerCalls = open(cancerCallsFile, 'r')
    normalCalls = open(normalCallsFile, 'r')
    germCalls   = open(germCallsFile, 'r')
    otherCalls  = open(otherCallsFile, 'r')
    refGenome   = pysam.Fastafile(refGenomeFile)

    cancerBreaksOut = open(args.outDirName + "/" + args.outBaseName + "/cancerbreaks.tab.txt", 'w') 
    normalBreaksOut = open(args.outDirName + "/" + args.outBaseName + "/normalbreaks.tab.txt", 'w')
    germBreaksOut   = open(args.outDirName + "/" + args.outBaseName + "/germlinebreaks.tab.txt", 'w')
    otherBreaksOut  = open(args.outDirName + "/" + args.outBaseName + "/uncategorizedbreaks.tab.txt", 'w')

    callSetListNames   = ('cancer', 'normal', 'germ','other')
    callSetListInFiles = (cancerCalls, normalCalls, germCalls, otherCalls)
    callSetListOutFiles = (cancerBreaksOut, normalBreaksOut, germBreaksOut, otherBreaksOut)

    for i in range(len(callSetListNames)):
        for line in callSetListInFiles[i]:
            col    = line.strip().split("\t")
            chr    = col[0]
            start  = int(col[1])
            end    = int(col[2])
            gene   = col[7]

            cancerCluster = fetchRegion(cancerBam,refGenome,int(args.maxReadLen),chr,start,end,gene,args.zeroChar,int(args.minClipQual),args.usechr)
            cancerCluster.type='CANCER'

            normalCluster = fetchRegion(normalBam,refGenome,int(args.maxReadLen),chr,start,end,gene,args.zeroChar,int(args.minClipQual),args.usechr)
            normalCluster.type='NORMAL'

            mergeCluster = mergeClusters(cancerCluster,normalCluster,txs)
            clusterout = mergeCluster.outstring()
            infodumpout = mergeCluster.infodump()

            callSetListOutFiles[i].write(line.strip("\n") + "\t" + clusterout + "\n" + infodumpout + "\n")
        callSetListInFiles[i].close()
        callSetListOutFiles[i].close()

if __name__ == '__main__':
    # commandline args
    parser = argparse.ArgumentParser(description='parse the output of discordant.py')
    parser.add_argument('-c', '--config', dest='configFileName', default='config.txt',
                        help='config file left by discordant.py')
    parser.add_argument('-o', '--outbasename', dest='outBaseName', required=True,
                        help='basename for output files')
    parser.add_argument('-d', '--outdirname', dest='outDirName', default='output',
                        help='output directory')
    parser.add_argument('-e-', '--eltfile', dest='eltFile', default='sumEltList.txt',
                        help='list of element families to include')
    parser.add_argument('-l', '--maxReadLen', dest='maxReadLen', default=100,
                        help='max read length in basepairs (default  100 bp)')
    parser.add_argument('-z', '--zerochar', dest='zeroChar', default='#',
                        help='for fastq quality scores, the character corresponding to zero (default #)')
    parser.add_argument('-g', '--refgenome', dest='refGenomeFile', required=True,
                        help='ref genome in fasta format, indexed with samtools faidx')
    parser.add_argument('-q', '--minclipqual', dest='minClipQual', default=30,
                        help='minimum avg. quality cutoff for trimmed region (default 30)')
    parser.add_argument('-m', '--mrnafile', dest='mrnaFastaFile', required=True,
                        help='directory of FASTA files with TE reference sequences in them, plus a config.txt file with ref names')
    parser.add_argument('--usechr', action="store_true", default=False,
                        help='set if reference genome uses "chr" prefix (default=False)')
    args = parser.parse_args()

    main(args)
