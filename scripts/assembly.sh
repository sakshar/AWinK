#!/bin/bash
# ------- PRE-REQUISITES
# Install the following packages in your conda environment before running this script
# hifiasm (channel: bioconda)
# command to install them:
# conda install -c bioconda hifiasm

READ_FILE="path to AWinK selected read file (FASTA)"
OUT_DIR="path to hifiasm output directory"
PREFIX="prefix for hifiasm outputs"
THREADS="provide the number of threads; for example 20"

hifiasm -o $OUT_DIR/$PREFIX  -t $THREADS $READ_FILE
awk '/^S/{print ">"$2;print $3}' $OUT_DIR/$PREFIX.bp.p_ctg.gfa > $OUT_DIR/$PREFIX.asm.fasta