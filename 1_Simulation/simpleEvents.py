# -*- coding: utf-8 -*-

import collections
import HTSeq
from sets import Set

SEgenes = open("SEgenes.txt")
RIgenes = open("RIgenes.txt")
A3genes = open("A3genes.txt")
A5genes = open("A5genes.txt")
otherevents = open("otherevents.txt")

SElist = list()
for line in SEgenes:
	SElist.append(line.rstrip())

RIlist = list()
for line in RIgenes:
	RIlist.append(line.rstrip())

A3list = list()
for line in A3genes:
	        A3list.append(line.rstrip())

A5list = list()
for line in A5genes:
        A5list.append(line.rstrip())

otherlist = list()
for line in otherevents:
        otherlist.append(line.rstrip())

versusSE = sum([RIlist, A3list, A5list, otherlist],[])
versusRI = sum([SElist, A3list, A5list, otherlist],[])
versusA3 = sum([RIlist, SElist, A5list, otherlist],[])
versusA5 = sum([RIlist, A3list, SElist, otherlist],[])

for item in versusSE:
        while item in SElist:
                SElist.remove(item)

for item in versusRI:
        while item in RIlist:
                RIlist.remove(item)

for item in versusA3:
        while item in A3list:
                A3list.remove(item)

for item in versusA5:
        while item in A5list:
                A5list.remove(item)

duplicates = list([x for x in SElist if SElist.count(x) > 1])

for item in duplicates:
        while item in SElist:
                SElist.remove(item)

duplicates = list([x for x in RIlist if RIlist.count(x) > 1])

for item in duplicates:
        while item in RIlist:
                RIlist.remove(item)

duplicates = list([x for x in A3list if A3list.count(x) > 1])

for item in duplicates:
        while item in A3list:
                A3list.remove(item)

duplicates = list([x for x in A5list if A5list.count(x) > 1])

for item in duplicates:
        while item in A5list:
                A5list.remove(item)

GENE = Set(["gene", "pseudogene", "transposable_element_gene"])
EXON = Set(["exon", "pseudogenic_exon"])
name_gff = "/scratch/leuven/313/vsc31305/IBP/TAIR10_GFF3_genes.gff"

num_lines = sum(1 for line in open(name_gff))

file_gff=open(name_gff,'r')
gff_file=HTSeq.GFF_Reader(file_gff)

count = 0
transcript = set()
lines = 0
gene_list=list()
for feature in gff_file:
        lines += 1
        if feature.type in GENE or lines == num_lines:
	        if len(transcript) == 2:
                        count += 1
                        gene_list.append(gene_cand.attr["ID"])
                gene_cand = feature
                transcript.clear()
        if feature.type in EXON:
                transcript.add(feature.attr["Parent"])


simpleEvents = open('simpleEvents.txt', 'w')

for item in SElist:
	if item in gene_list:
		print >> simpleEvents, item + "\t" + "SE" + "\t" + "2"		

for item in RIlist:
        if item in gene_list:
                print >> simpleEvents, item + "\t" + "RI" + "\t" + "2"

for item in A3list:
        if item in gene_list:
		print >> simpleEvents, item + "\t" + "A3" + "\t" + "2"

for item in A5list:
        if item	in gene_list:
		print >> simpleEvents, item + "\t" + "A5" + "\t" + "2"                                                                                                                                                                

simpleEvents.close()
