#! /bin/bash

# Execute this script from the the root/top-level project directory.
#
# The script requires a single parameter, which is the name of
# the database to be created and loaded. You can provide a path
# such as /database/synteny.db in which case the database file will
# be created in the /database directory.

DATA_DIR='data/'
SRC_DIR='src/'

# Loads:
# -genes and exons for Human, Mouse, and Rat
echo Loading genes and exons data
python ${SRC_DIR}import_genes_exons.py $1 ${DATA_DIR}Human_GenomeFeatures_Synteny.gff3 9606 -c
python ${SRC_DIR}import_genes_exons.py $1 ${DATA_DIR}Mouse_GenomeFeatures_Synteny.gff3 10090
python ${SRC_DIR}import_genes_exons.py $1 ${DATA_DIR}Rat_GenomeFeatures_Synteny.gff3 10116

# - syntenic blocks for Human, Mouse, and Rat
echo Loading syntenic blocks data for Human and Rat
python ${SRC_DIR}import_synteny_blocks.py $1 ${DATA_DIR}HumanRat_Blocks_Synteny.txt -c
echo Loading syntenic blocks data for Mouse and Human
python ${SRC_DIR}import_synteny_blocks.py $1 ${DATA_DIR}MouseHuman_Blocks_Synteny.txt
echo Loading syntenic blocks data for Rat and Mouse
python ${SRC_DIR}import_synteny_blocks.py $1 ${DATA_DIR}RatMouse_Blocks_Synteny.txt

# - cytobands for Human, Mouse, and Rat
echo Loading cytogenentic bands loci data
python ${SRC_DIR}import_cytogenetic_loc.py $1 ${DATA_DIR}Human_CytoBand_Synteny.gff3 9606 -c
python ${SRC_DIR}import_cytogenetic_loc.py $1 ${DATA_DIR}Mouse_CytoBand_Synteny.gff3 10090
python ${SRC_DIR}import_cytogenetic_loc.py $1 ${DATA_DIR}Rat_CytoBand_Synteny.gff3 10116

# - QTL data for Mouse and Rat
echo Loading QTL data
python ${SRC_DIR}import_features.py $1 ${DATA_DIR}Mouse_QTL_Synteny.gff3 10090 -c
python ${SRC_DIR}import_features.py $1 ${DATA_DIR}Rat_QTL_Synteny.gff3.gz 10116

# - GWAS variants for Human
echo Loading GWAS variants data
python ${SRC_DIR}import_variants.py $1 ${DATA_DIR}Human_GWAS_Synteny.vcf 9606 -c

# Loads: GO, MP, DO
echo Loading ontology data
python ${SRC_DIR}import_ontology.py $1

# Load the homologs (syntenic regions) from a file.
echo Loading homolog data for Mouse and Human
python ${SRC_DIR}import_homologs.py $1 ${DATA_DIR}MouseHuman_Homologs_Synteny.tsv -c
echo Loading homolog data for Rat and Human
python ${SRC_DIR}import_homologs.py $1 ${DATA_DIR}RatHuman_Homologs_Synteny.tsv
echo Loading homolog data for Rat and Mouse
python ${SRC_DIR}import_homologs.py $1 ${DATA_DIR}RatMouse_Homologs_Synteny.tsv

echo FINISHED