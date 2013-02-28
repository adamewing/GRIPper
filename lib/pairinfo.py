#!/usr/bin/env python

"""
 Copyright (C) 2012 by Adam Ewing (adam.ewing@gmail.com)

 Released under the MIT license, see LICENSE.txt
"""

import math

class Pair:
    """Stores information for a genome/element read pair"""
    def __init__(self,chrom,readPos,readStrand,mateElt,mateChrom,mateStrand,
                 matePos,eltStart,eltEnd,eltStrand,eltFullLen,source,peakIndex,geneName):
        self.chrom      = chrom
        self.readPos    = readPos
        self.readStrand = readStrand
        self.mateElt    = "_".join(mateElt.split(' '))
        self.mateChrom  = mateChrom
        self.mateStrand = mateStrand
        self.matePos    = matePos
        self.eltStart   = eltStart
        self.eltEnd     = eltEnd
        self.eltStrand  = eltStrand
        self.eltFullLen = eltFullLen # length of gene in this case
        self.source     = source # keep track of genome
        self.peakIndex  = peakIndex
        self.geneName   = "_".join(geneName.split(' '))

    def relStrand(self):
        """return the strand of the ref element relative to the alignment"""
        if self.eltStrand == self.mateStrand:
            return '+'
        else:
            return '-'

    def relPos(self):
        """return the position of the mate read relative to the ref element coords"""
        refLen   = self.eltEnd - self.eltStart
        refLeft  = self.matePos - self.eltStart
        refRight = self.eltEnd - self.matePos
        if self.eltStrand == '+':
            if (int(self.eltFullLen) - refRight) < 0:
                return 0
            else:
                return int(self.eltFullLen) - refRight
        else: # - strand
            if (int(self.eltFullLen) - refLeft) < 0:
                return 0
            else:
                return int(self.eltFullLen) - refLeft

class Peak:
    """Stores objects of class Pair and provides functions for analysis"""
    def __init__(self):
        self.pairs = []
        self.npairs = 0 
    def addpair(self,pair):
        self.pairs.append(pair)
        self.npairs = self.npairs + 1
    def __iter__(self):
        return self
    def next(self):
        if self.npairs == 0:
            self.index = len(self.pairs)
            raise StopIteration
        self.npairs = self.npairs - 1
        return self.pairs[self.npairs]

    def getBreakpoint(self):
        """find breakpoint (where ref strands switch)"""
        plusScore  = [0] * len(self.pairs) 
        minusScore = [0] * len(self.pairs)
        index = 0
        for pair in self.pairs:
            if pair.readStrand == '+':
                if index == 0:
                    plusScore[0]  = 1
                    minusScore[0] = 0
                else:
                    plusScore[index]  = plusScore[index-1] + 1
                    minusScore[index] = minusScore[index-1] - 1
                if minusScore[index] < 0:
                    minusScore[index] = 0

            if pair.readStrand == '-':
                if index == 0:
                    plusScore[0]  = 0
                    minusScore[0] = 1
                else:
                    minusScore[index] = minusScore[index-1] + 1
                    plusScore[index]  = plusScore[index-1] - 1
                if plusScore[index] < 0:
                    plusScore[index] = 0

            index = index + 1

        leftHiScorePlus  = 0
        leftHiScoreMinus = 0
        leftHiIndexPlus  = 0
        leftHiIndexMinus = 0
        
#       print "forward result (going left to right across the peak):"
        for i in range(0,index):
            if plusScore[i] > leftHiScorePlus:
                leftHiScorePlus = plusScore[i]
                leftHiIndexPlus = i
            if minusScore[i] > leftHiScoreMinus:
                leftHiScoreMinus = minusScore[i]
                leftHiIndexMinus = i
            
#       print "%d %d %d" % (i, plusScore[i], minusScore[i])
#       print "leftHiScorePlus: %d" % leftHiScorePlus + " index: %d" % leftHiIndexPlus
#       print "leftHiScoreMinus: %d" % leftHiScoreMinus + " index: %d" % leftHiIndexMinus

        # reset and run in the other direction across the peak
        plusScore  = [0] * len(self.pairs)
        minusScore = [0] * len(self.pairs)
        index = len(self.pairs) - 1
        for pair in reversed(self.pairs):
            if pair.readStrand == '+':
                if index == len(self.pairs) - 1:
                    plusScore[index]  = 1
                    minusScore[index] = 0
                else:
                    plusScore[index]  = plusScore[index+1] + 1
                    minusScore[index] = minusScore[index+1] - 1
                if minusScore[index] < 0:
                    minusScore[index] = 0

            if pair.readStrand == '-':
                if index == len(self.pairs) - 1:
                    plusScore[index]  = 0
                    minusScore[index] = 1
                else:
                    minusScore[index] = minusScore[index+1] + 1
                    plusScore[index] = plusScore[index+1] - 1
                if plusScore[index] < 0:
                    plusScore[index] = 0

            index = index - 1

        rightHiScorePlus  = 0
        rightHiScoreMinus = 0
        rightHiIndexPlus  = 0
        rightHiIndexMinus = 0

#        print "reverse result (going right to left across the peak):"
        for i in range(0,len(self.pairs)):
            if plusScore[i] > rightHiScorePlus:
                rightHiScorePlus = plusScore[i]
                rightHiIndexPlus = i
            if minusScore[i] > rightHiScoreMinus:
                rightHiScoreMinus = minusScore[i]
                rightHiIndexMinus = i

#        print "%d %d %d" % (i, plusScore[i], minusScore[i])
#        print "rightHiScorePlus: %d" % rightHiScorePlus + " index: %d" % rightHiIndexPlus
#        print "rightHiScoreMinus: %d" % rightHiScoreMinus + " index: %d" % rightHiIndexMinus

        if (leftHiIndexPlus-rightHiIndexMinus < leftHiIndexMinus-rightHiIndexPlus):
            myBreakPoint = BreakPoint(leftHiIndexPlus,
                                      rightHiIndexMinus,
                                      self.pairs[leftHiIndexPlus].readPos,
                                      self.pairs[rightHiIndexMinus].readPos,
                                      '+', '-',self.npairs,self)
            return myBreakPoint
	else: # will need to filter out bad breakpoints elsewhere since they're all reported
            myBreakPoint = BreakPoint(leftHiIndexMinus,
                                      rightHiIndexPlus,
                                      self.pairs[leftHiIndexMinus].readPos,
                                      self.pairs[rightHiIndexPlus].readPos,
                                      '-','+',self.npairs,self)
            return myBreakPoint

    def countWrongStrand(self,indexStart,indexEnd,expectedStrand):
        """
        Return the number of reads that are not on the expected strand between indexStart
        and indexEnd, inclusive.
        """
        wsCount = 0
        for index in range(indexStart,indexEnd+1):
            if (self.pairs[index].readStrand != expectedStrand):
                wsCount = wsCount + 1
        return wsCount

    def countWrongMateStrand(self,indexStart,indexEnd,expectedStrand):
        wsCount = 0
        for index in range(indexStart,indexEnd+1):
            if (self.pairs[index].relStrand() != expectedStrand):
                wsCount = wsCount + 1
        return wsCount

    def meanDist(self,indexStart,indexEnd):
        """
        Compute the mean distance between reads across a range indices. Note that since a range
        statement range(indexStart,indexEnd+1) is used, indexStart and indexEnd will both be 
        included in the mean. The distance for each index is the difference in position between
        the read with that index and the read with index-1 (so it doesn't make sense to call this
        with index 0)
        """
        distSum = 0
        for index in range(indexStart,indexEnd+1):
            distSum = distSum + (self.pairs[index].readPos - self.pairs[index-1].readPos)
        if (indexEnd + 1 - indexStart) < 1:
            return 0
        else:
            return distSum/(indexEnd + 1 - indexStart)

    def meanMatePos(self,indexStart,indexEnd):
        """
        Get mean position of reads within the insertion for a range of indices (inclusive)
        needs the expected length of a full-length insertion to determine the relative position
        inside the insertion from the genomic coords
        """
        posSum = 0
        for index in range(indexStart,indexEnd+1):
            posSum = posSum + self.pairs[index].matePos
        if (indexEnd + 1 - indexStart) < 1:
            return 0
        else:
            return posSum/(indexEnd + 1 - indexStart)

    def maxMatePos(self,indexStart,indexEnd):
        """as meanMatePos, but returns maximum position in element"""
        maxPos = self.pairs[indexStart].matePos
        for index in range(indexStart,indexEnd+1):
            myPos = self.pairs[index].matePos
            if (myPos > maxPos):
                maxPos = myPos
        return maxPos

    def minMatePos(self,indexStart,indexEnd):
        """as meanMatePos, but returns minimum position in element"""
        minPos = self.pairs[indexStart].matePos
        for index in range(indexStart,indexEnd+1):
            myPos = self.pairs[index].matePos
            if (myPos < minPos):
                minPos = myPos
        return minPos

    def majorityExonType(self,indexStart,indexEnd):
        """Determine whether most exons are first (F), last (L), internal (I) or unitary (O)"""
        exoncount = {}
        for index in range(indexStart,indexEnd+1):
            gex  = self.pairs[index].mateElt.split('.')
            gene = gex[0]
            exon = gex[len(gex)-1]
            if exon != 'L' and exon != 'F' and exon !='O':
                exon = 'I'
            if exoncount.has_key(exon):
                exoncount[exon] += 1
            else:
                exoncount[exon] = 1
        majority = ''
        max = 0
        for exon,count in exoncount.iteritems():
            if count > max:
                max = count
                majority = exon
        return majority

    def hasPairs(self):
        """Return true if this peak has anything in it (sanity check)."""
        if len(self.pairs) < 1:
            return 0
        else:
            return 1

    def getChr(self):
        """return chromosome name"""
        if (self.pairs[0]):
            return self.pairs[0].chrom
        else:
            return 0

    def getSources(self):
        """
        # returns a dictionary of the sources of reads for this peak and the number
        # of reads each source was responsible for
        """
        sourceDict = dict()
        for pair in self.pairs:
            if sourceDict.has_key(pair.source):
                sourceDict[pair.source] = sourceDict[pair.source] + 1
            else :
                sourceDict[pair.source] = 1

        return sourceDict

    def primaryMateElt(self):
        """Returns the element family to which the majority of mates are mapped"""
        eltDict = {}
        for pair in self.pairs:
            if eltDict.has_key(pair.mateElt):
                eltDict[pair.mateElt] += 1
            else:
                eltDict[pair.mateElt] = 1
        maxEltCount = 0
        maxEltName = ''
        for elt,count in eltDict.iteritems():
            if count > maxEltCount:
                maxEltCount = count
                maxEltName  = elt
        return maxEltName

    def eltStrand(self, indexLeft, indexRight):
        """Infers the element strand given the indices of the reads flanking the breakpoint"""
        leftStrand  = self.pairs[indexLeft].mateStrand
        rightStrand = self.pairs[indexRight].mateStrand
        if leftStrand == '-' and rightStrand == '+':
            return '+'
        if leftStrand == '+' and rightStrand == '-':
            return '-'
        if leftStrand == rightStrand:
            if leftStrand == '+': # inversion due to twin priming (Ostertag et al.)
                return 'i'
            if leftStrand == '-': # ??
                return 'x'

    def getIndex(self):
        return int(self.pairs[0].peakIndex)

    def numExons(self):
        """
        pseudogene-specific
        max number of exon types in __one gene__, an exon is considered present if there
        are at least two reads representing it
        types are: F=first, L=last, O=only, I=internal (converted from an exon #)
        """
        genecount = {}
        seen      = {}
        for pair in self.pairs:
            gex  = pair.mateElt.split('.')
            gene = gex[0]
            exon = gex[len(gex)-1]
            if exon != 'L' and exon != 'F' and exon !='O':
                exon = 'I'

            exname = gene + "." + exon

            if seen.has_key(exname):
                seen[exname] += 1
            else:
                seen[exname] = 1 

        for exname,count in seen.iteritems():
            if count > 1:
                (gene,exon) = exname.split('.')
                if genecount.has_key(gene):
                    genecount[gene] += 1
                else:
                    genecount[gene] = 1

        max = 0
        for gene,count in genecount.iteritems():
            if count > max:
                max = count
        return max

    def exonTypes(self):
        """Returns a dictionary of exon types present"""
        types = {} 
        for pair in self.pairs:
            exonparse = pair.mateElt.split('.')
            types[exonparse[len(exonparse)-1]] = 1
        return types.keys()

class BreakPoint:
    """Stores information about a breakpoint consisting of a 5' side and a 3' side"""
    def __init__(self,leftindex,rightindex,leftpos,rightpos,leftstrand,rightstrand,npairs,peak):
        self.leftindex   = leftindex
        self.rightindex  = rightindex
        self.leftpos     = leftpos
        self.rightpos    = rightpos
        self.leftstrand  = leftstrand
        self.rightstrand = rightstrand
        self.npairs      = npairs
        self.peak        = peak

    def leftMeanDist(self):
        """Function to find the mean distance between reads on the left side of the breakpoint."""
        return self.peak.meanDist(1,self.leftindex)
    def rightMeanDist(self):
        """Function to find the mean distance between reads on the right side of the breakpoint."""
        return self.peak.meanDist(self.rightindex+1,self.npairs-1)
    
    def leftWrongStrand(self):
        """How many reads in left peak at the insertion site are on the wrong strand"""
        return self.peak.countWrongStrand(0,self.leftindex,self.leftstrand)

    def rightWrongStrand(self):
        """How many reads in the right peak at the insertion site are on the wrong strand"""
        return self.peak.countWrongStrand(self.rightindex,self.npairs-1,self.rightstrand)

    def leftMateWrongStrand(self):
        """How many reads in left peak at the parent gene are on the wrong strand"""
        leftMateStrand = self.peak.pairs[self.leftindex].relStrand()
        return self.peak.countWrongMateStrand(0,self.leftindex,leftMateStrand)

    def rightMateWrongStrand(self):
        """How many reads in the right peak at the parent gene are on the wrong strand"""
        rightMateStrand = self.peak.pairs[self.rightindex].relStrand()
        return self.peak.countWrongMateStrand(self.rightindex,self.npairs-1,rightMateStrand)

    def minPeakPos(self):
        """leftmost peak coordinatei at insertion site"""
        return self.peak.pairs[0].readPos
    def maxPeakPos(self):
        """rightmost peak coordinate at insertion site"""
        return self.peak.pairs[self.npairs-1].readPos
    def leftMateMeanPos(self):
        """mean position of the left reads on the parent gene"""
        return self.peak.meanMatePos(0,self.leftindex)
    def rightMateMeanPos(self):
        """mean position of the right reads on the parent gene"""
        return self.peak.meanMatePos(self.rightindex,self.npairs-1)
    
    def leftMajorExonType(self):
        """return majority exon type from left read cluster, exon types are F,I,O,L as described elsewhere"""
        return self.peak.majorityExonType(0,self.leftindex)
    def rightMajorExonType(self):
        """return majority exon type from right read cluster, exon types are F,I,O,L as described elsewhere"""
        return self.peak.majorityExonType(self.rightindex,self.npairs-1)

    # min positions of reads in elements
    def leftMateMinPos(self):
        """minimum (leftmost) position of reads mapped to parent gene in left read cluster"""
        return self.peak.minMatePos(0,self.leftindex)
    def rightMateMinPos(self):
        """minimum (leftmost) position of reads mapped to parent gene in right read cluster"""
        return self.peak.minMatePos(self.rightindex,self.npairs-1)

    # max positions of reads in elements
    def leftMateMaxPos(self):
        """maximum (rightmost) position of reads mapped to parent gene in left read cluster"""
        return self.peak.maxMatePos(0,self.leftindex)
    def rightMateMaxPos(self):
        """maximum (rightmost) position of reads mapped to parent gene in right read cluster"""
        return self.peak.maxMatePos(self.rightindex,self.npairs-1)

    def eltStrandFromStrand(self):
        """element strand from read strands (can be 'i' if inverted or 'x' if FUBAR)"""
        return self.peak.eltStrand(self.leftindex,self.rightindex)

    def eltStrandFromReads(self):
        """gets element strand from where the reads are mapped within the element"""
        minside='left'
        maxside='right'
        if (self.rightMateMinPos() < self.leftMateMinPos()):
            minside='right'
        if (self.leftMateMaxPos() > self.rightMateMaxPos()):
            maxside='left'

        if (minside == maxside):
            return 'x'
        else:
            if (minside=='left' and maxside=='right'):
                return '+'
            if (minside=='right' and maxside=='left'):
                return '-'

    def sourceList(self):
        return self.peak.getSources().keys()

    def mateOverlapsInsert(self):
        """
        Some putative pseudogene insertion sites are very close to the mapping position of their mates.
        these probably represent deletions, inversions, or some other confounding SV and not an actual
        pseudogene insertion
        """
        window = 5000 # arbitrary
        for pair in self.peak.pairs:
            if pair.mateChrom == pair.chrom: # pair.chrom should be the same for all pairs in the cluster...
                if pair.matePos > self.minPeakPos()-window and pair.matePos < self.maxPeakPos()+window:
                    return True
        return False

    def bestStrand(self):
        """
        uses the other two strand calling function to make a best guess as to the orientation
        of the element, returns '.' by default
        """
        sFS = self.eltStrandFromStrand()
        sFR = self.eltStrandFromReads()
        if sFS == sFR and sFS != 'x':
            return sFS
        if (sFS == 'i' or sFS == 'x') and sFR != 'x':
            return sFR
        if sFR == 'x' and sFS != 'x' and sFS != 'i':
            return sFS
        if sFR != 'x' and sFS !='x' and sFS != 'i' and sFR != sFS:
            return sFS
        return '.'

    def eltMinPos(self):
        """minimum position within the insertion (parent gene)"""
        return min(self.leftMateMinPos(),self.rightMateMinPos())
    def eltMaxPos(self):
        """maximum position within the insertion (parent gene)"""
        return max(self.leftMateMaxPos(),self.rightMateMaxPos())

    def judgeExonComposition(self):
        """
        Processed pseudogenes should have at least a 3' UTR, and the left and right aligned regions
        (either side of the breakpoint) should have a different majority exon type. Since some genes
        have very short 5' exons and a longer 5' or 3' UTR, we can also accept mappings that have a
        wide enough spacing between the mapped mates minimum and maximum positions relative to the
        mappings to the insertion site.
        """
        leftMajor = self.leftMajorExonType()
        rightMajor = self.rightMajorExonType()
        if leftMajor == rightMajor: # left and right should be in differnet exons - at least mostly
            mateLeftPos  = self.peak.pairs[self.leftindex].matePos
            mateRightPos = self.peak.pairs[self.rightindex].matePos
            if abs(mateLeftPos-mateRightPos) > 2*abs(self.rightpos-self.leftpos):
                #print "abs(%d-%d) > abs(%d,%d)" % (mateLeftPos, mateRightPos, self.rightpos, self.leftpos)
                return True
            return False
        # somebody has to be the 3' end
        if leftMajor != 'L' and rightMajor != 'L':
            return False
        return True

    def isValid(self):
        """
        returns whether a breakpoint meets the following conditions:
        genomic 'anchor' alignments are either +/- or -/+ (5'/3')
        reads are within one position of each other
        peak size cutoff
        """
        if math.fabs(self.leftindex-self.rightindex) > 1:
#            print "inv rule 1"
            return False
        if self.leftstrand != '+' or self.rightstrand != '-':
#            print "inv rule 2"
            return False
        # require two reads on each junction
        if self.leftindex < 1 or self.npairs-self.rightindex-1 < 1:
#            print "inv rule 3"
            return False
        # requirements for mean inter-read distances inside of a peak
        if self.leftMeanDist() > 200 or self.rightMeanDist() > 200:
#            print "inv rule 4"
            return False
        # less than 10% of reads on either side can have an inconsistent strand
#        if ((float(self.leftWrongStrand()) / float(self.leftindex+1) > .25) or 
#            (float(self.rightWrongStrand()) / float(self.rightindex-self.leftindex+1) > .25)):
#            print "inv rule 5"
#            return False
        if float(self.leftWrongStrand() + self.rightWrongStrand())/float(self.npairs) > .1:
#            print "inv rule 6"
            return False
        if self.mateOverlapsInsert():
#            print "inv rule 9"
            return False
        if self.rightpos-self.leftpos > 300:
#            print "inv rule 10"
            return False
        if not self.judgeExonComposition():
#            print "inv rule 11"
            return False
        return True




