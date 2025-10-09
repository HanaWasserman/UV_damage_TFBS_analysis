import pandas as pd
import os
from Bio.Seq import Seq
import logomaker as lm
import numpy as np

def prep_motif(tf_path):
    """
    Make TFBS logo
    """
    seqs = pd.read_table(tf_path, header=None)
    sqs = seqs[0].str.upper()
    sqs = sqs.loc[~sqs.str.contains('N')]
    rev_complements = sqs.apply(lambda x: str(Seq(x).complement()))
    counts_mat = lm.alignment_to_matrix(sqs.tolist())
    counts_mat.index = counts_mat.index - (20 / 2) + 1
    rev_counts_mat = lm.alignment_to_matrix(rev_complements.tolist())
    rev_counts_mat.index = counts_mat.index  # Align the indices with original
    #Normalize to probability
    num_sequences = len(sqs)
    pwm = counts_mat.apply(lambda x: x / num_sequences, axis=1)
    rev_pwm = rev_counts_mat.apply(lambda x: x / num_sequences, axis=1)
    ic_matrix = lm.transform_matrix(pwm, from_type='probability', to_type='information')
    rev_ic_matrix = lm.transform_matrix(rev_pwm, from_type='probability', to_type='information')
    #Multiply reverse complement matrix by -1 for the negative y-axis
    rev_ic_matrix *= -1
    return ic_matrix, rev_ic_matrix

def get_tf_motif_seqs(tf, tf_file, tf_window_size, tf_dir, genome_dir):
    slop_size = tf_window_size-10
    os.system(f'bedtools slop -i {tf_dir}/{tf_file}.bed -g {genome_dir}/hg19/hg19.fa.fai -l -{slop_size} -r -{slop_size} > {tf_dir}/{tf}_top_20.bed')
    os.system(f'bedtools getfasta -fi {genome_dir}/hg19/hg19.fa -bed {tf_dir}/{tf}_top_20.bed -fo {tf_dir}/{tf}_top_20_seqs -s -bedOut')

def curate_pvals(df, strand, window_size):
    df_strand = df.loc[(df['pos'].between(-window_size,window_size))&(df['strand']==strand)].copy()
    df_strand['qval'] = df_strand['qval'].replace(0,1e-10)
    return df_strand

def merge_reps(df1, df2, qv_cutoff=.05, cpd_cutoff=30, zscore="zscore_scal"):
    df1_df2 = pd.merge(df1, df2[['pos', 'enriched','strand','TF','pval','qval',zscore]], on=['TF','pos','strand'], suffixes=['_cell', '_nDNA'])
    df1_df2['FDR'] = False
    df1_df2.loc[(df1_df2['qval_cell']<qv_cutoff)&(df1_df2['count']>=cpd_cutoff)&(((df1_df2['enriched_cell']==df1_df2['enriched_nDNA'])&(df1_df2['qval_nDNA']>=0))|(df1_df2['enriched_cell']!=df1_df2['enriched_nDNA'])), 'FDR'] = True
    df1_df2['qval_log'] = float(0)
    df1_df2.loc[df1_df2['enriched_cell']==False,'qval_log'] = np.log10(df1_df2['qval_cell'])
    df1_df2.loc[df1_df2['enriched_cell']==True,'qval_log'] = np.log10(df1_df2['qval_cell'])*-1
    df1_df2['zscore_diff'] = df1_df2[f'{zscore}_cell']-df1_df2[f'{zscore}_nDNA']
    return df1_df2

def prep_reps(cell_pvals, nDNA_pvals, window_size, qv_cutoff=.05, cpd_cutoff=30):
    cell_pvals_same = curate_pvals(cell_pvals, 'same', window_size)
    cell_pvals_opp = curate_pvals(cell_pvals, 'opp', window_size)
    nDNA_pvals_same = curate_pvals(nDNA_pvals, 'same', window_size)
    nDNA_pvals_opp = curate_pvals(nDNA_pvals, 'opp', window_size)
    return merge_reps(cell_pvals_same, nDNA_pvals_same, qv_cutoff, cpd_cutoff), merge_reps(cell_pvals_opp, nDNA_pvals_opp, qv_cutoff, cpd_cutoff)
