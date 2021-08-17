# -*- coding: utf-8 -*-

"""

Beagle requires atleast one marker per contig (which the program calls windows).
This quick script generates a list of markers for Beagle to exclude.
The output file of this script can direclty be used as input for Beagle.

The output is of the format LOCUS:POSTION, one marker per line.

No CLI for this script yet - just change the file names in the variable(s) below 
and run. 
"""

from collections import Counter

VCF_file = 'raw.g5mac3dp3.recode.vcf' 
OUTPUT_file = 'excludemarkers.txt'


#step1: get all the locuses

locus = []
with open(VCF_file, 'r') as in_f:
    line = in_f.readline()
    while line:
        if line.startswith('#') != True:
            l = line.split('\t')[0]
            locus.append(l)
        line = in_f.readline()


# create a dictionary using counting the number of entries per locus
counts = Counter(locus)
singles = [l for l in counts.keys() if counts[l] == 1] #separate loci occuring only once
pos = []

# get position co-ordinates for singles
with open(VCF_file, 'r') as in_f:
    line = in_f.readline()
    while line:
        if line.startwith('#') != True:
            l = line.split('\t')[0]
            if l in singles:
                pos.append(line.split('\t')[1])
        line = in_f.readline()

# write to output file.
with open(OUTPUT_file, 'w') as out_f:
    for i in range(len(singles)):
        out_f.write(singles[i]+":"+pos[i]+"\n")