#!/bin/env python

import logging, os, re, argparse, peakparser, pickreads, mergepairs

# +1 if a > b
# -1 if a < b
#  0 if a == b
def cmpChrPosList(a,b):

    chrA = a[0].lstrip('chr')
    # handle X/Y
    chrA = re.sub("X","2000",chrA)
    chrA = re.sub("Y","3000",chrA)
    # handle 2a/2b (chimp)
    chrA = re.sub("2a","4000",chrA)
    chrA = re.sub("2b","5000",chrA)
    # handle '_random' chromosomes
    if re.search("_random", chrA):
        chrA = int(chrA.rstrip("_random")) + 1000
    if chrA < ':' and chrA > '0': chrA = int(chrA)

    chrB = b[0].lstrip('chr')
    chrB = re.sub("X","2000",chrB)
    chrB = re.sub("Y","3000",chrB)
    chrB = re.sub("2a","4000",chrB)
    chrB = re.sub("2b","5000",chrB)
    if re.search("_random", chrB):
        chrB = int(chrB.rstrip("_random")) + 1000
    if chrB < ':' and chrB > '0': chrB = int(chrB)

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

# sorts files whose first two columns are chromosome and position, respectively.
def sortChrPos(file):
    f = open(file,'r')
    bedlist = []
    for bedline in f:
        bedlist.append(bedline.rstrip().split("\t"))
    return sorted(bedlist, cmp=cmpChrPosList)

# concatenates files, sorts, and outputs to outfile
def mergeChrPosFiles(infiles,outfile,maxDist):
    bedlist = []
    for file in infiles:
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
        outlist = mergepairs.reindexReads(sorted(bedlist, cmp=cmpChrPosList), maxDist)
    else:
        outlist = sorted(bedlist, cmp=cmpChrPosList)

    f = open(outfile, 'w')
    if re.search("track", bedhead):
        f.write("track\tname=%s\tvisibility=2\titemRgb=On db=hg19\n" % os.path.basename(outfile))
    for outline in outlist: #sorted(bedlist, cmp=cmpChrPosList):
        outstring = "\t".join(outline)
        f.write(outstring + "\n")
    f.close()

def main(args):

    # create output directory
    pickreads.prepOutDir(args.outBaseName,args.outDirName,args.overwrite)

    pickreads.checkfile(args.sampleListFile)

    sampleList = open(args.sampleListFile, 'r')

    readFileNames = []
    bamFileNames  = []
    insertSizes   = []
    readLengths   = []
    lastConfig    = None

    for sampleLine in sampleList:
        if not re.search("^#", sampleLine):
            (sampleBam,sampleSubDir,refGenome) = sampleLine.strip().split()
            peakparser.checkOutDir(sampleSubDir,args.outDirName)
            configPath = args.outDirName + "/" + sampleSubDir + "/" + args.configFileName

            pickreads.checkfile(configPath)
            configDict = peakparser.readConfig(configPath,sampleSubDir,args.outDirName)
            lastConfig = configDict

            insertSizes.append(int(configDict['insertSize']))
            readLengths.append(int(configDict['readLength']))

            readFileName = args.outDirName + "/" + sampleSubDir + "/" + sampleSubDir + ".readpairs.txt"
            readFileNames.append(readFileName)
            bamFileNames.append(configDict['bamFileName'])

    maxDist = max(insertSizes) + 2*max(readLengths)        

    # merge readfiles
    outReadFileName = args.outDirName + "/" + args.outBaseName + "/" + args.outBaseName + ".readpairs.txt"
    mergeChrPosFiles(readFileNames,outReadFileName,maxDist)

    # write new config file
    configPath = args.outDirName + "/" + args.outBaseName + "/" + args.configFileName
    configDict = lastConfig

    bfnum = 0
    for bamFileName in bamFileNames:
        bfvname = "bamFileName" + str(bfnum)
        configDict[bfvname] = bamFileName
        bfnum += 1

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
    parser.add_argument('-o', '--outbasename', dest='outBaseName', required=True,
                        help='basename for merged output')
    parser.add_argument('-d', '--outdirname', dest='outDirName', default='output',
                        help='output directory')
    parser.add_argument('-c', '--config', dest='configFileName', default='config.txt',
                        help='filname for config file left by pickreads.py')
    parser.add_argument('-s', '--sampleList', dest='sampleListFile', required=True)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()

    main(args)


