#!/bin/bash

#Paths
GENOME_DIR="/Users/hanawasserman/PycharmProjects/genomes"
MAIN_DIR="/Users/hanawasserman/PycharmProjects/UV_damage_TFBS_analysis"
DATA_DIR="/Users/hanawasserman/PycharmProjects/CPDdata"

#Parameters
TF_WINDOW_SIZE=100
KMER_LENGTH=6
DAMAGE_FORMATION_WINDOW=30
DAMAGE_REPAIR_WINDOW=30
BOOTSTRAP_SAMPS=100
TFBS_BUFFER=5
SIM_MODE="cumulative"

#Files
ARCHETYPE_FILE="archetype_motifs_intersect_WT_CSB_hg19_idr_conservative_summits_150bp_intergenic.bed"
DHS_FILE="WT_CSB_hg19_idr_conservative_summits_150bp_intergenic_sorted_merged.bed"
MOTIF_CLUSTER_FILE="TF_cluster_motifs_sample.csv"
EXPERIMENTS=("WT_CSB_6J_0hr_CPD_1bp_sorted" "WT_CSB_6J_6hr_CPD_1bp_sorted" "XPC_12J_NakedDNA_S22_CPD_1bp_sorted")
TIMECOURSE_EXP1="WT_CSB_6J_0hr_CPD_1bp_sorted"
TIMECOURSE_EXP2="WT_CSB_6J_6hr_CPD_1bp_sorted"