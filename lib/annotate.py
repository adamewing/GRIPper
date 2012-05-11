#!/bin/env python

"""
 Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)

 Released under the MIT license, see LICENSE.txt
"""

import pysam
import sys
import re

def uniqann(annList):
    uniq = {}
    for ann in annList:
        uniq[ann] = 1
    return uniq.keys()

def checkfile(fname):
    try:
        open(fname)
    except IOError as e:
        print "can't find file: " + fname
        sys.exit()


class annotator:
    def __init__(self,tabixfile,name,genome):
        checkfile(tabixfile)
        
        self.tabixfn = tabixfile  # location of tabix file
        self.tabix = pysam.Tabixfile(tabixfile, 'r')

        self.name    = name   # name for annotation
        self.genome  = genome # genome assembly

        self.tabixContigs = {}
        for contig in self.tabix.contigs:
            self.tabixContigs[contig] = 1

    # returns a list of annotations
    def annotate(self,bedline,genome):
        c = bedline.rstrip().rsplit("\t")
        chr   = c[0]
        start = c[1]
        end   = c[2]

        if not re.search('chr',chr):
            raise LookupError("chromosome names must start with chr: " + chr)
            return []
        if (self.genome != genome):
            raise LookupError("tried to compare a %s bedfile to a %s annotation." % (genome,self.genome))
            return []
        else:
            annotations = []
            if (chr and start and end):
                if self.tabixContigs.has_key(chr):
                    tabixTupleParse = self.tabix.fetch(reference=chr,
                                                       start=int(start),
                                                       end=int(end),
                                                       parser=pysam.asTuple())
                    for tabixTuple in tabixTupleParse:
                        annotations.append(tabixTuple[3])
                    return uniqann(annotations)
                else:
                    return []
            else:
                raise LookupError("can't find chr,start,end. File must be tab-delimited")
                return []

if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit()
    a = annotator(sys.argv[1],'default','hg19')
    f = open(sys.argv[2], 'r')
    for line in f:
        if not re.search("^track", line):
            for annotation in a.annotate(line,'hg19'):
                print annotation
    f.close()
