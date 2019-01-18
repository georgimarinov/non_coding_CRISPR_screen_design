# Non-coding CRISPR-Cas9 screen design tools

These scripts can be used to design and analyze non-coding CRISPR screens with GuideScan specificity scores. 

**Library design:** The `extractGuidesFromGuideScan.py` script will take a set of regions and output the desired number of sgRNAs within those regions, based on either the number of off-targets (with 2-3 mismatches) or the GuideScan specificity scores.

**Screen analysis:** The `GuidesPerRegionFromWholeGenomeGuideScan.py` script can be used to obtain specificity scores for a set of sgRNAs based on their spacer sequences (provided those guides have been included in the GuideScan database). This can be used for pre-existing libraries that were not initially generated with GuideScan.

**Input files:** sgRNA scores were downloaded from the GuideScan webtool (www.guidescan.com) and can be found here (for hg38 and mm10):

http://mitra.stanford.edu/kundaje/marinovg/non_coding_CRISPR_screen_design/
