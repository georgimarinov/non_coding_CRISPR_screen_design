# non_coding_CRISPR_screen_design

These scripts can be used to design non-coding CRIPSR screens based on GuideScan specificity scores. 

The extractGuidesFromGuideScan.py script will take a set of regions and output the desired number of sgRNAs within those regions, based on either the number of off-targets or the GuideScan CFD sgRNA specificity scores.

The GuidesPerRegionFromWholeGenomeGuideScan.py script can be used to obtain CFD scores for a set of sgRNAs (provided those guides have been included in the GuideScan database). 

The GuideScan files used as input can be found here (for hg20 and mm10):

http://mitra.stanford.edu/kundaje/marinovg/non_coding_CRISPR_screen_design/
