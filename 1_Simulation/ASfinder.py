# -*- coding: utf-8 -*-

'''
Program to find AS genes in a genome annotation file (gff). Output is a bed file containing the regions of annotated AS genes. 
Author: Carmen Bravo, Carlos de Lannoy, Emma Houben
Version: 1.0
'''	

# Load the hpc modules before running this:
#module load HTSeq
#module load IPython/3.0.0-foss-2014a-Python-2.7.6
	

import HTSeq
from sets import Set

GENE = Set(["gene", "transposable_element_gene","pseudogene"])
EXON = Set(["exon", "pseudogenic_exon"])
name_gff = raw_input('Select gff file:')

num_lines = sum(1 for line in open(name_gff))

file_gff=open(name_gff,'r')
gff_file=HTSeq.GFF_Reader(file_gff)

count = 0
transcript = set()
gene_list = list()
chrom_list = list()
start_list = list()
end_list = list()
strand_list = list()
lines = 0
for feature in gff_file:
	lines += 1
	if feature.type in GENE:
		if len(transcript) > 1 or lines == num_lines:
			count += 1
			gene_list.append(gene_cand.attr["ID"])
			chrom_list.append(gene_cand.iv.chrom)
			start_list.append(gene_cand.iv.start)
			end_list.append(gene_cand.iv.end)
			strand_list.append(gene_cand.iv.strand)
			#print(gene_cand.__dict__)
		gene_cand = feature
		transcript.clear()
	if feature.type in EXON:
		transcript.add(feature.attr["Parent"])
print 'Number of genes', count
print 'Length of AS genes', len(gene_list)

ASids = open('ASgenome.txt', 'w')
i=0
while i < len(gene_list):
	print >> ASids, str(chrom_list[i]) + "\t" + str(start_list[i]) + "\t" +  str(end_list[i]) + "\t" + gene_list[i] + "\t" + "0" + "\t" + str(strand_list[i])
	i += 1

ASregs = open('ASgenome_reg.bed', 'w')
i=0
while i < len(gene_list):
        print >> ASregs, str(chrom_list[i]) + "\t" + str(start_list[i]) + "\t" +  str(end_list[i]) + "\t" + gene_list[i] + "\t" + "0" + "\t" + str(strand_list[i])
        i += 1
file_gff.close()
ASids.close()
ASregs.close()
