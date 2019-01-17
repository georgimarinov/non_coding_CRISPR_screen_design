##################################
#                                #
# Last modified 2018/08/28       #
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
import os
from sets import Set
import time

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s guidescan_input wanted_guides wanted_guides_fieldID [-splitby string] [-minGuideLength N] [-5pG]' % sys.argv[0]
        print '\tthe script will print to stdout by default'
        print '\tNote: use - for stdin for the guidescan input; the script will read compressed files for the wanted guides automatically'
        print '\tthe script will look for guides 20bp long, lower the -minGuideLength parameter for shorter guides'
        print '\tuse the [-5pG] option to enable looking for mismatched 5pG'
        print '\tAssumed GuideScan format:'
        print '\t\tchromosome,target site start coordinate,target site end coordinate,gRNA,cutting efficiency score,cutting specificity score,strand,offtargets sum,offtargets summary,annotation,gRNA label'
        sys.exit(1)

    GS = sys.argv[1]
    wanted = sys.argv[2]
    fieldID = int(sys.argv[3])

    do5G = False
    if '-5pG' in sys.argv:
        do5G = True
#        print 'will look for mismatched 5pG bases'

    minGL = 20
    if '-minGuideLength' in sys.argv:
        minGL = int(sys.argv[sys.argv.index('-minGuideLength') + 1])
#        print 'will look for guides as short as', minGL

    SplitString = '\t'
    if '-splitby' in sys.argv:
        SplitString = sys.argv[sys.argv.index('-splitby') + 1]

    WantedGuideDict = {}

    if wanted.endswith('.bz2'):
        cmd = 'bzip2 -cd ' + wanted
    elif wanted.endswith('.gz'):
        cmd = 'gunzip -c ' + wanted
    elif wanted.endswith('.zip'):
        cmd = 'unzip -p ' + wanted
    else:
        cmd = 'cat ' + wanted
    p = os.popen(cmd, "r")
    line = 'line'
    while line != '':
        line = p.readline().strip()
        fields = line.split('\t')
        if line == '':
            break
        if line.startswith('#'):
            continue
        fields = line.strip().split(SplitString)
        sgRNA = fields[fieldID]
        WantedGuideDict[sgRNA] = 1

    if GS == '-':
        linelist = sys.stdin
    else:
        linelist = open(GS)
    for line in linelist:
        if line.startswith('#'):
            continue
        fields = line.strip().split(',')
        if len(fields) < 3:
            continue
        sgRNA = fields[3]
        if WantedGuideDict.has_key(sgRNA):
            print line.strip()
        else:
            Found = False
            while not Found and len(sgRNA) >= minGL:
                sgRNA = sgRNA[0:-1]
                if WantedGuideDict.has_key(sgRNA):
                    Found = True
                    print line.strip()
                if not Found and do5G:
                    if WantedGuideDict.has_key('G' + sgRNA):
                        Found = True
                        print line.strip()
                    elif WantedGuideDict.has_key('G' + sgRNA[1:]):
                        Found = True
                        print line.strip()

run()
