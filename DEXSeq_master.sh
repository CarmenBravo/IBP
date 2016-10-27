#!/bin/bash

# DEXseq master script
#
# Runs two pythonscripts to consecutively "flatten" the gtf-file and convert bam-files into count files.
# Then runs Rscript to perform analysis and parse results into a table of genes and p-values.

DATA_PATH=$1 #Specify where to look for necessary data; BAM-files, gtf/gff-file
NCORES=$2 # Specify number of cores to use in R-part of script
PY_PATH=$(Rscript DEXSeq_findPythonScripts.R | grep -Po '(?<=\").*(?=\")') # Retrieve location of python files
GTF=$(find ${DATA_PATH}/feature_files -name *.gtf)
NbGTF=$(find ${DATA_PATH}/feature_files -name *.gtf| wc -l)

# Check if only one GTF-file exists in the file structure. If not: throw error
if [ $NbGTF != "1" ]; then
echo "Error: feature_files folder must contain 1 gtf file"
exit 1
fi

echo $GTF
# Make flattened GFF-file out of GTF file.
python ${PY_PATH}/dexseq_prepare_annotation.py ${DATA_PATH}/feature_files/TAIR10_GTF_genes.gtf ${DATA_PATH}/feature_files/DEXSeq_FlattenedFeatureFile.gff

for file in ${DATA_PATH}/bam_files/*.bam; do
python ${PY_PATH}/dexseq_count.py -f bam ${DATA_PATH}/feature_files/DEXSeq_FlattenedFeatureFile.gff ${file}  ${DATA_PATH}/DEXSeq_output/HTSeqCount_files/${file}_count.txt
done

Rscript DEXSeqScript.R $DATA_PATH $NCORES
