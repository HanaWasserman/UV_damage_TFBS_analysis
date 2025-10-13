import sys
import os
import pandas as pd
pd.options.mode.chained_assignment = None
from utils.pipeline_functions import rev_complement, get_tf
import progressbar

def prep_tfbs(tf_cluster, tf_len, tf, tf_path, genome_path, archetype_file, tf_window_size, archetype_path):
    """
    Subset binding site windows for TF of interest
    """
    if not os.path.exists(tf_path):
        os.system(f"mkdir {tf_path}")

    window_shft_sz = tf_window_size - (tf_len // 2)
    gen_fa_fai, gen_fa = f'{genome_path}/hg19/hg19.fa.fai', f'{genome_path}/hg19/hg19.fa'
    os.system(f"grep -w {tf_cluster} {archetype_path}/{archetype_file} > {tf_path}/{tf}_tmp_unsorted.bed")
    os.system(f"bedtools sort -i {tf_path}/{tf}_tmp_unsorted.bed > {tf_path}/{tf}_tmp_sorted.bed")
    tf_bed_unfilt = pd.read_table(f'{tf_path}/{tf}_tmp_sorted.bed', header=None, usecols=[0, 1, 2, 3, 4, 5], names=['chrom', 'start', 'end', 'misc', 'score', 'strand'])
    tf_bed_unfilt['wind_sz'] = tf_bed_unfilt['end']-tf_bed_unfilt['start']
    tf_bed_filt = tf_bed_unfilt.loc[tf_bed_unfilt['wind_sz']==tf_len].copy()
    tf_bed_filt[['chrom', 'start', 'end', 'misc', 'score', 'strand']].to_csv(f"{tf_path}/{tf}.bed", header=None, index=False, sep='\t')
    os.remove(f"{tf_path}/{tf}_tmp_unsorted.bed")
    os.remove(f"{tf_path}/{tf}_tmp_sorted.bed")

    def slop_it(window_shft_sz):
        """
        Create temporary .bed file with modified for coordinates using 'slop' for correct fasta retrieval
        """
        os.system(f"bedtools slop -i {tf_path}/{tf}.bed -g {gen_fa_fai} -l {window_shft_sz} -r {window_shft_sz} > {tf_path}/{tf}_{window_shft_sz}_tmp_slop.bed")
        """Retrieve fasta sequences for coordinates in modified .bed file"""
        os.system(f"bedtools getfasta -fi {gen_fa} -bed {tf_path}/{tf}_{window_shft_sz}_tmp_slop.bed -fo {tf_path}/{tf}_{window_shft_sz}_tmp_slop_seq -s -bedOut")

        tf_bed = pd.read_table(f'{tf_path}/{tf}_{window_shft_sz}_tmp_slop.bed', header=None, usecols=[0,1,2,4,5], names=['chrom','start','end','score','strand'])
        seqs = pd.read_table(f"{tf_path}/{tf}_{window_shft_sz}_tmp_slop_seq", header=None)
        tf_bed['seq'] = seqs[0].str.upper()
        tf_bed_complete = rev_complement(tf_bed)
        os.remove(f"{tf_path}/{tf}_{window_shft_sz}_tmp_slop.bed")
        os.remove(f"{tf_path}/{tf}_{window_shft_sz}_tmp_slop_seq")
        tf_bed_reordered = tf_bed_complete[['chrom','start','end','seq','seq_rv','strand', 'score']].copy()
        tf_bed_reordered.to_csv(f"{tf_path}/{tf}_{tf_window_size}_seq.bed", header=None, index=False, sep='\t')

        tf_score_0 = tf_bed_reordered['score'].quantile(0 / 2, interpolation='nearest')
        tf_score_1 = tf_bed_reordered['score'].quantile(1 / 2, interpolation='nearest')
        tf_score_2 = tf_bed_reordered['score'].quantile(2 / 2, interpolation='nearest')
        tf_1 = tf_bed_reordered.loc[tf_bed_reordered['score'].between(tf_score_0, tf_score_1)].copy().reset_index(drop=True)
        tf_2 = tf_bed_reordered.loc[tf_bed_reordered['score'].between(tf_score_1, tf_score_2)].copy().reset_index(drop=True)
        tf_1.to_csv(f"{tf_path}/{tf}_bottom_{tf_window_size}_seq.bed", header=None, index=False, sep='\t')
        tf_2.sort_values(by='score', ascending=False)[:10000].to_csv(f"{tf_path}/{tf}_top_{tf_window_size}_seq_tmp.bed", header=None, index=False, sep='\t')
        os.system(f"bedtools sort -i {tf_path}/{tf}_top_{tf_window_size}_seq_tmp.bed > {tf_path}/{tf}_top_{tf_window_size}_seq.bed")
        os.remove(f"{tf_path}/{tf}_top_{tf_window_size}_seq_tmp.bed")
        #tf_2.to_csv(f"{tf_path}/{tf}_top_{tf_window_size}_seq.bed", header=None, index=False, sep='\t')

    slop_it(window_shft_sz)

def roi_kmers(tf_file, kmer, tf_path, offset=2):
    '''
    Count kmers in aggregated TF binding site window positions by NYYN sequence
    '''
    yys = ['TT', 'CT', 'TC', 'CC']
    sites, window_sz = get_tf(f'{tf_path}/{tf_file}_seq', ['chr', 'start', 'end', 'seq', 'rv_seq', 'strand'])
    def count_kmers(windows, plus_strand):
        """
        Boolean:param plus_strand: True means TFBS is on the + strand and False means TFBS is on - strand
        """
        kmer_list = []
        for seq in progressbar.progressbar(windows):
            seq_l = len(seq)
            for i in range(seq_l - (kmer-1)):
                seq_n = seq[i:i + kmer]
                if 'N' in seq_n:
                    continue
                if seq_n[offset:(offset+2)] in yys:
                    if plus_strand:
                        kmer_list.append({'pos':i, 'seq':seq_n})
                    else:
                        kmer_list.append({'pos':seq_l - i - kmer, 'seq':seq_n})
        kmer_df = pd.DataFrame(kmer_list)
        kmer_df['count'] = 1
        kmer_agg = kmer_df.groupby(by = ['pos', 'seq'], as_index=False).sum()
        return kmer_agg

    print(f'Counting {tf_file} {kmer}mers...')
    p_kmers = count_kmers(sites['seq'], True)
    m_kmers = count_kmers(sites['rv_seq'], False)
    p_kmers.to_csv(f"{tf_path}/{tf_file}_{kmer}mers_ROI_plus.csv", index=False)
    m_kmers.to_csv(f"{tf_path}/{tf_file}_{kmer}mers_ROI_minus.csv", index=False)

if __name__ == "__main__":
    main_dir = sys.argv[1]
    genome_path = sys.argv[2]
    archetype_file = sys.argv[3]
    tf_window_size = int(sys.argv[4])
    kmer = int(sys.argv[5])
    tf_cluster = sys.argv[6]
    tf_len = int(sys.argv[7])
    tf = sys.argv[8]
    tf_path = f'{main_dir}/results/TFBS'
    data_path = f'{main_dir}/data'

    print(f"Processing {tf_cluster}...")
    prep_tfbs(tf_cluster, tf_len, tf, tf_path, genome_path, archetype_file, tf_window_size, data_path)
    roi_kmers(f"{tf}_top_{tf_window_size}", kmer, tf_path)
    roi_kmers(f"{tf}_bottom_{tf_window_size}", kmer, tf_path)