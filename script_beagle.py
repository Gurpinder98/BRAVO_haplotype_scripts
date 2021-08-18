# -*- coding: utf-8 -*-

import numpy as np
from collections import OrderedDict

#TODO make a CLI
GFF_FILE = "Brassica_napus.AST_PRJEB5043_v1.44.sorted.gff3"
#GENE_TO_LOOK = ['BnaC02g00490D'] #temporary arragement
GENE_TO_LOOK = ['BnaA02g00370D']
VCF_FILE = 'test2.vcf'

# create a structure {LKXXXX: {Gene1:(start,stop), Gene2:(start,strop)}, LKXXX:{{}}.. }

Genes = {}
genes_found = 0
with open(GFF_FILE, "r") as in_f:
    line = in_f.readline()
    while line:
        line = line.rstrip('\n')
        if line.startswith("#") != True:
            contig = line.split('\t')[0]
            locus_type = line.split('\t')[2]
            if locus_type == 'gene': #Can also look at UTRs, exons etc. 
                if contig not in Genes.keys():
                    Genes[contig] = {}
                gene_name = line.split("\t")[8].split(";")[1].lstrip('Name=')
                if gene_name in GENE_TO_LOOK: #only if gene is in GENE_TO_LOOK it is added
                    (start_pos, stop_pos) = (line.split("\t")[3], line.split("\t")[4])
                    Genes[contig][gene_name] = (start_pos, stop_pos)
                    genes_found +=1
        line = in_f.readline()


#clearing all contigs with no genes
all_contigs = list(Genes.keys())        
for contig in all_contigs:
    if Genes[contig] == {}:
        Genes.pop(contig)

print("Successfully read {}.\n Found {} contigs for {}/{} genes.".format(GFF_FILE, len(Genes.keys()), genes_found, len(GENE_TO_LOOK)))

# complete_allele_patterns: {SLXXXX:{LKXXXX:{'0001','10101'}, LKXX:..}, SLXX..}
complete_allele_patterns = {}
# super_contigs: { LKXXXX: OrderedDict{POS1: (ref, alt), POS2: (ref, alt)..}, LKXXX..}
#order of posttion values in dict is crucial for this - so LKXXX dicts are all ordered.  
super_contigs = {}

with open(VCF_FILE, 'r') as in_f:
    line = in_f.readline()
    
    while line:
        line = line.rstrip('\n')
        if line.startswith('#CHROM'):
            samples = line.split('\t')[9:]
            print("Reading the {} file, {} samples found.".format(VCF_FILE, len(samples)))
            for indv in samples:
                complete_allele_patterns[indv] = {}
        
        if line.startswith('#') != True:
            locus = line.split('\t')[0]
            if locus in list(Genes.keys()): #only if locus is of interest it is considered.
                super_contigs[locus] = OrderedDict()
                for indv in complete_allele_patterns.keys():
                    complete_allele_patterns[indv][locus] = ['','']
        line = in_f.readline()

print("All datastructures constructed, found {} contigs from VCF file for specified genes.".format(len(super_contigs.keys())))


with open(VCF_FILE, 'r') as in_f:
    line = in_f.readline()
    while line:
        line = line.rstrip('\n')
        if line.startswith('#') != True:
            line_array = line.split('\t')
            locus = line_array[0]
            position = line_array[1]
            ref = line_array[3]
            alt = line_array[4]
            #variant_type = line_array[7].split(';')[-1].lstrip('TYPE=')
            if locus in Genes.keys():
                super_contigs[locus][position] = (ref, alt)
                
                #for genotype
                all_samples_genotypes = [data for data in line_array[9:]]
                assert len(all_samples_genotypes) == len(complete_allele_patterns.keys())
                for i in range(len(all_samples_genotypes)):
                    complete_allele_patterns[list(complete_allele_patterns.keys())[i]][locus][0] = complete_allele_patterns[list(complete_allele_patterns.keys())[i]][locus][0] + all_samples_genotypes[i].split("|")[0]
                    complete_allele_patterns[list(complete_allele_patterns.keys())[i]][locus][1] = complete_allele_patterns[list(complete_allele_patterns.keys())[i]][locus][1] + all_samples_genotypes[i].split("|")[1]
                
        line = in_f.readline()
        
print("All data loaded from {}".format(VCF_FILE))        


def alleles_slicer(start_stop, allele_seqs, positions):
        
    """
        Takes whole contig allele sequences and outputs sliced sequences.
        
        input:
            start_stop: a tuple of ints 
            allele_seqs: a tuple of two alleles patterns. a array/tuple of strings.
            positions: a list of positions for every allele in haplotype sequence.
        
        returns:
            a tupple of sliced allele seqs
    """
    
    try:
        assert len(allele_seqs[0]) == len(positions)
        assert len(allele_seqs[0]) == len(allele_seqs[1])
    except:
        AssertionError
        print("Lengths of position array and both allele strings do not match.")
    pos_array = []
    sliced_allele_0 = ''
    sliced_allele_1 = ''
    for i in range(len(positions)):
        if int(positions[i]) >= int(start_stop[0]) and int(positions[i]) <= int(start_stop[1]):
            pos_array.append(positions[i])
            sliced_allele_0 = sliced_allele_0 + allele_seqs[0][i]    
            sliced_allele_1 = sliced_allele_1 + allele_seqs[1][i]
    
    return (pos_array, [sliced_allele_0, sliced_allele_1])
    

# final data structure {Gene:{sample:[marker_arr, alllele1, allele2]}}
final_data = {}
for super_locus in super_contigs.keys():
    super_locus_dict = {}
    position_array = list(super_contigs[super_locus].keys())
    for gene in Genes[super_locus].keys():
        current_gene_dict = {}
        start_stop = Genes[super_locus][gene]
        for sample in complete_allele_patterns.keys():
            current_gene_dict[sample] = {} #fill in the innermost dict first
            
            raw_allele_seqs = complete_allele_patterns[sample][super_locus]
            pos_array, sliced_allele_seqs = alleles_slicer(start_stop, raw_allele_seqs, position_array)
            
            #marker format : Locus:POS:Ref allele:alt_allele1,alt_allele2
            marker_array = [super_locus+':'+pos+':'+super_contigs[super_locus][pos][0]+':'+super_contigs[super_locus][pos][1] for pos in pos_array]
            current_gene_dict[sample] = [marker_array, sliced_allele_seqs[0], sliced_allele_seqs[1]]
        final_data[gene] = current_gene_dict
        
       
        
#TODO Distance Matrix

samples = list(final_data[GENE_TO_LOOK[0]].keys())
print("{} samples in samples list".format(len(samples)))
Similarity_Matrix = np.zeros((len(samples), len(samples)))

def distance_calculator(seqA_array, seqB_array):
    """
    Return distance between two loci based on two set of allele patterns

    Parameters
    ----------
    seqA_array : List of two strings
        eg. ['101001', '101101']
    seqB_array : List of two strings
        eg. ['101001', '101101']
        

    Returns
    -------
    float
        distance between two samples, averaged. 
        Scoring scheme: +1 if both allele patterns differ
                        +0.5 if only one is different
    """
    try:
        assert len(seqA_array[0]) == len(seqB_array[0])
        assert len(seqA_array[1]) == len(seqB_array[1])
    except:
        AssertionError
        print("Sequence length mismatch.")
    pattern0 = zip(seqA_array[0], seqB_array[0])
    difference0 = sum([1 for allele in pattern0 if allele[0] != allele[1]])
    pattern1 = zip(seqA_array[1], seqB_array[1])
    difference1 = sum([1 for allele in pattern1 if allele[0] != allele[1]])

    return (difference0 + difference1)/2

for vertical in range(Similarity_Matrix.shape[0]):
    for horizontal in range(Similarity_Matrix.shape[1]):
        Similarity_Matrix[vertical][horizontal] = distance_calculator(final_data[GENE_TO_LOOK[0]][samples[vertical]][1:],final_data[GENE_TO_LOOK[0]][samples[horizontal]][1:])

# create a dendrogram
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

import matplotlib.pyplot as plt

with open("sample_names.txt", "r") as in_f:
    lines = in_f.readlines()

sample_names_dict = {}
for line in lines[1:]:
    sample_names_dict[line.split("\t")[0]] = line.split("\t")[1]
sample_line_names = [sample_names_dict[s] for s in samples]

dists = squareform(Similarity_Matrix)
linkage_matrix = linkage(dists, "single")
plt.figure(figsize=(15,7))
plt.box(False)
dendrogram(linkage_matrix, labels=sample_line_names)
plt.title(GENE_TO_LOOK[0])
plt.savefig(GENE_TO_LOOK[0]+".pdf", dpi=300)




#TODO final output
