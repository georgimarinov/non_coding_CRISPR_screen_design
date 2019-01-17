##################################
#                                #
# Last modified 2019/01/16       #
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
import random
import os
import gzip
from sets import Set
import time

def run():

    if len(sys.argv) < 7:
        print 'usage: python %s guidescan_prefix region chrFieldID leftFieldID rightFieldID N_guides_per_region outfile [-addAdapter left right] [-simpleOverlap] [-extendRegion left right] [-flankControls distance radius GTF exon|CDS] [-sortByCFD]' % sys.argv[0]
        print '\tguidescan prefix example:'
        print '\t\t\thg38-male'
        print 'where files are named as follows:'
        print '\t\t\thg38-male.chr16.csv.gz'
        print '\tBy default the script will pick guides only among the set of guides with the fewest predicted off-target sites for each region'
        print '\tuse the [-sortByCFD] option to sort by the CFD score'
        print '\tBy default the script will only ouput Cas9 guides that cut within the region; if you want all guides overlapping the region, use the [-simpleOverlap] option'
        print '\tUse the [-flankControls] to pick control guides nearby within a distance of the a minimum distance provided by the distance and radius parameters and not overlapping either the exons or the CDSs in the specified GTF file'
        sys.exit(1)

    GS = sys.argv[1]
    regions = sys.argv[2]
    chrFieldID = int(sys.argv[3])
    leftFieldID = int(sys.argv[4])
    rightFieldID = int(sys.argv[5])
    N = int(sys.argv[6])
    outfilename = sys.argv[7]

    doCFDsort = False
    if '-sortByCFD' in sys.argv:
        doCFDsort = True

    doER = False
    if '-extendRegion' in sys.argv:
        doER = True
        ERLeft = int(sys.argv[sys.argv.index('-extendRegion') + 1])
        ERRight = int(sys.argv[sys.argv.index('-extendRegion') + 2])
        print 'will extend regions by', ERLeft, 'bp on the left side, and by', ERRight, 'on the right'

    doAddAd = False
    if '-addAdapter' in sys.argv:
        doAddAd = True
        AdapterLeft = sys.argv[sys.argv.index('-addAdapter') + 1]
        AdapterRight = sys.argv[sys.argv.index('-addAdapter') + 2]
        print 'will add adapter sequences to guides:', AdapterLeft, AdapterRight	

    doSO = False
    if '-simpleOverlap' in sys.argv:
        doSO = True
        print 'will output all guides overlapping the region'

    doFlanks = False
    if '-flankControls' in sys.argv:
        doFlanks = True
        D = int(sys.argv[sys.argv.index('-flankControls') + 1])
        R = int(sys.argv[sys.argv.index('-flankControls') + 2])
        GTF = sys.argv[sys.argv.index('-flankControls') + 3]
        ExonOrCDS = sys.argv[sys.argv.index('-flankControls') + 4]
        print 'will also pick one guide at a distance of ~', D, 'while excluding guides possibly targetting', ExonOrCDS, 'sequences'

    RegionDict = {}

    if regions.endswith('.gz'):
        linelist = gzip.open(regions)
    else:
        linelist = open(regions)
    for line in linelist:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chr = fields[chrFieldID]
        left = int(fields[leftFieldID])
        right = int(fields[rightFieldID])
        if doER:
            left = left - ERLeft
            right = right + ERRight
        if RegionDict.has_key(chr):
            pass
        else:
            RegionDict[chr] = {}
        for i in range(left,right):
            RegionDict[chr][i] = 1
        if doFlanks:
            for i in range(left - D - R,left - D + R):
                RegionDict[chr][i] = 1
            for i in range(right + D - R,right + D + R):
                RegionDict[chr][i] = 1

    if doFlanks:
        ExonDict = {}
        j=0
        lineslist = open(GTF)
        for line in lineslist:
            j+=1
            if j % 100000 == 0:
                print j, 'lines processed'
            if line.startswith('#'):
                continue
            fields=line.strip().split('\t')
            if fields[2] != ExonOrCDS:
                continue
            chr = fields[0]
            if RegionDict.has_key(chr):
                pass
            else:
                continue
            left=int(fields[3])
            right=int(fields[4])
            orientation=fields[6]
            if ExonDict.has_key(chr):
                pass
            else:
                ExonDict[chr] = {}
            for i in range(left,right):
                ExonDict[chr][i] = 1

    print 'finished inputting regions'

    for chr in RegionDict.keys():
        input = GS + '.' + chr + '.csv.gz'
        print chr, GS + '.' + chr + '.csv.gz'
        try:
            linelist = open(input)
        except:
            print 'could not find file', input, 'skipping'
            continue
        if input.endswith('.bz2'):
            cmd = 'bzip2 -cd ' + input
        elif input.endswith('.gz'):
            cmd = 'gunzip -c ' + input
        elif input.endswith('.zip'):
            cmd = 'unzip -p ' + input
        else:
            cmd = 'cat ' + input
        p = os.popen(cmd, "r")
        line = 'line'
        while line != '':
            line = p.readline().strip()
            if line == '':
                break
            if line.startswith('chromosome,target site'):
                continue
            if ',' not in line:
                continue
            fields = line.strip().split(',')
            C = fields[0]
            if C != chr:
                print 'mismatch between file name and chromosome content, exiting'
                print C, chr, GS + '.' + chr + '.csv.gz'
                sys.exit(1)
            left = int(fields[1])
            right = int(fields[2])
            strand =  fields[6]
            gRNA = fields[3]
            cutting_efficiency_score = fields[4]
            cutting_specificity_score = fields[5]
            offtargets_sum = int(fields[7])
            offtargets_summary = fields[8]
            if doCFDsort:
                guide = (cutting_specificity_score,offtargets_sum,offtargets_summary,chr,left,right,gRNA,cutting_efficiency_score,strand)
            else:
                guide = (offtargets_sum,offtargets_summary,chr,left,right,gRNA,cutting_efficiency_score,cutting_specificity_score,strand)
            if doSO:
                if strand == '+':
                    cut = right - 3
                if strand == '-':
                    cut = left + 3
                if RegionDict[chr].has_key(right):
                    if RegionDict[chr][right] == 1:
                        RegionDict[chr][right] = []
                    RegionDict[chr][right].append(guide)
                elif RegionDict[chr].has_key(left):
                    if RegionDict[chr][left] == 1:
                        RegionDict[chr][left] = []
                    RegionDict[chr][left].append(guide)
                else:
                    continue
            else:
                if strand == '+':
                    cut = right - 3
                if strand == '-':
                    cut = left + 3
                if RegionDict[chr].has_key(cut):
                    pass
                else:
                    continue
                if RegionDict[chr][cut] == 1:
                    RegionDict[chr][cut] = []
                RegionDict[chr][cut].append(guide)

    outfile = open(outfilename, 'w')

    SeenDict = {}

    if regions.endswith('.gz'):
        linelist = gzip.open(regions)
    else:
        linelist = open(regions)
    for line in linelist:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chr = fields[chrFieldID]
        left = int(fields[leftFieldID])
        right = int(fields[rightFieldID])
        if doER:
            left = left - ERLeft
            right = right + ERRight
        guides = []
        for i in range(left,right):
            if RegionDict[chr][i] != 1:
                for g in RegionDict[chr][i]:
                    guides.append(g)
        guides = list(Set(guides))
        guides.sort()
        if doCFDsort:
            guides.reverse()
        if len(guides) >= 1:
            pass
        else:
            continue
        for i in range(min(len(guides),N)):
            if doCFDsort:
                (cutting_specificity_score,offtargets_sum,offtargets_summary,chr,gleft,gright,gRNA,cutting_efficiency_score,strand) = guides[i]
            else:
                (offtargets_sum,offtargets_summary,chr,gleft,gright,gRNA,cutting_efficiency_score,cutting_specificity_score,strand) = guides[i]
            if SeenDict.has_key(gRNA):
                pass
            else:
                outline = chr + ':' + str(left) + '-' + str(right) + '__'
                outline = outline + chr + ':' + str(gleft) + '-' + str(gright) + ':' + strand + ',' + str(cutting_efficiency_score) + ',' + str(cutting_specificity_score) + ',' + str(offtargets_sum) + ';' + offtargets_summary
                if doAddAd:
                    gRNA = AdapterLeft + gRNA + AdapterRight
                outline = outline + ',' + gRNA
                outfile.write(outline + '\n')
                SeenDict[gRNA] = 1
        if doFlanks:
            guides = []
            for i in range(left - D - R,left - D + R):
                if RegionDict[chr][i] != 1:
                    if ExonDict.has_key(chr):
                        if ExonDict[chr].has_key(i):
                            continue
                if RegionDict[chr][i] != 1:
                    for g in RegionDict[chr][i]:
                        guides.append(g)
            guides = list(Set(guides))
            guides.sort()
            if doCFDsort:
                guides.reverse()
            if len(guides) >= 1:
                if doCFDsort:
                    (cutting_specificity_score,offtargets_sum,offtargets_summary,chr,gleft,gright,gRNA,cutting_efficiency_score,strand) = guides[0]
                else:
                    (offtargets_sum,offtargets_summary,chr,gleft,gright,gRNA,cutting_efficiency_score,cutting_specificity_score,strand) = guides[0]
                if SeenDict.has_key(gRNA):
                    pass
                else:
                    outline = chr + ':' + str(left) + '-' + str(right) + '_leftFlank__'
                    outline = outline + chr + ':' + str(gleft) + '-' + str(gright) + ':' + strand + ',' + str(cutting_efficiency_score) + ',' + str(cutting_specificity_score) + ',' + str(offtargets_sum) + ',' + offtargets_summary
                    if doAddAd:
                        gRNA = AdapterLeft + gRNA + AdapterRight
                    outline = outline + ',' + gRNA
                    outfile.write(outline + '\n')
                    SeenDict[gRNA] = 1
            guides = []
            for i in range(right + D - R,right + D + R):
                if RegionDict[chr][i] != 1:
                    if ExonDict.has_key(chr):
                        if ExonDict[chr].has_key(i):
                            continue
                if RegionDict[chr][i] != 1:
                    for g in RegionDict[chr][i]:
                        guides.append(g)
            guides = list(Set(guides))
            guides.sort()
            if len(guides) >= 1:
                if doCFDsort:
                    (cutting_specificity_score,offtargets_sum,offtargets_summary,chr,gleft,gright,gRNA,cutting_efficiency_score,strand) = guides[0]
                else:
                    (offtargets_sum,offtargets_summary,chr,gleft,gright,gRNA,cutting_efficiency_score,cutting_specificity_score,strand) = guides[0]
                if SeenDict.has_key(gRNA):
                    pass
                else:
                    outline = chr + ':' + str(left) + '-' + str(right) + '_rightFlank__'
                    outline = outline + chr + ':' + str(gleft) + '-' + str(gright) + ':' + strand + ',' + str(cutting_efficiency_score) + ',' + str(cutting_specificity_score) + ',' + str(offtargets_sum) + ',' + offtargets_summary
                    if doAddAd:
                        gRNA = AdapterLeft + gRNA + AdapterRight
                    outline = outline + ',' + gRNA
                    outfile.write(outline + '\n')
                    SeenDict[gRNA] = 1

    outfile.close()

run()
