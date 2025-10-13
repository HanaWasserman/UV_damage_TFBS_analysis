import pandas as pd
import os
import sys
import progressbar
from utils.pipeline_functions import adjust_kmer, filter_kmer, get_tf

def proc_TFBS_dams(tf_file, kmer, exp, data_path, output_path, tf_path, genome_path, offset=1):

    if not os.path.exists(output_path):
        os.system(f"mkdir {output_path}")

    f_path = f"{output_path}/{tf_file}_intersect_{exp}_{kmer}mer"
    fasta = f"{genome_path}/hg19/hg19.fa"
    def prep_damage(strand_txt):
        print(f'Preparing: {f_path}')
        os.system(f"bedtools intersect -u -a {data_path}/{exp}_{strand_txt}_agg_filt.bed -b {tf_path}/{tf_file}.bed > {f_path}_{strand_txt}strand.bed")
        df = pd.read_table(f'{f_path}_{strand_txt}strand.bed', header=None, usecols=[0, 1, 2, 3], names=['chr', 'st', 'end', 'counts'])
        return df

    prep_df_plus = prep_damage('plus')
    prep_df_minus = prep_damage('minus')

    adj_df_plus = adjust_kmer(prep_df_plus, -(offset+1), (offset+2), f_path, 'plus', ' -s', fasta)
    adj_df_minus = adjust_kmer(prep_df_minus, -(offset+1), (offset+2), f_path, 'minus', '', fasta)

    print(f'Filtering: {exp}')
    filter_kmer(adj_df_plus, offset, 'plus', f'{f_path}')
    filter_kmer(adj_df_minus, offset, 'minus', f'{f_path}')

    os.remove(f"{f_path}_minusstrand.bed")
    os.remove(f"{f_path}_plusstrand.bed")
    os.remove(f"{f_path}_minus_tmp_tmp_seqs")
    os.remove(f"{f_path}_plus_tmp_tmp_seqs")
    os.remove(f"{f_path}_tmp_tmp.bed")

def count_damages(tf_file, kmer, exp, output_path, tf_path):

    f_txt = f"{output_path}/{tf_file}_intersect_{exp}_{kmer}mer_dipy_proc"
    chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    sites, window_sz = get_tf(f'{tf_path}/{tf_file}', ['chr', 'start', 'end', 'seq', 'rv_seq', 'strand'])

    def get_damages_from_file(fname):
        print("Loading " + fname)
        return pd.read_table(fname, header=None, usecols=[0, 1, 3, 4], names=['chr','start','counts','seq'])

    def count_it(chr, windows, wind_sz, damages_df, plus_strand):
        """
        Iterate through each TFBS window and count / record CPDs for each position in the window
        """
        damages_df_len = damages_df.shape[0]
        d_seqs = list()
        i = 0
        for bs_start, bs_end, bs_seq in windows.iloc:
            if bs_end-bs_start > wind_sz:
                continue
            while True:
                if i > damages_df_len-1:
                    break
                if damages_df[i]['start'] > bs_start:
                    if i > 0:
                        i = i - 1
                        continue
                    else:
                        break
                else:
                    break
            while i < damages_df_len:
                damage, count, seq = damages_df[i]
                if damage <= bs_start:
                    i = i + 1
                    continue
                if damage >= bs_end - 1:
                    break
                if plus_strand:
                    d_seqs.append((chr, damage, damage - bs_start, seq, count))
                elif not plus_strand:
                    d_seqs.append((chr, damage, bs_end - damage - kmer, seq, count))
                else:
                    print(f'CPD did not match on either strand: {damage} {seq}')
                i = i + 1
        return d_seqs

    p_damages = get_damages_from_file(f'{f_txt}_plus.bed')
    m_damages = get_damages_from_file(f'{f_txt}_minus.bed')

    """
    Initialize empty output array for total CPD counts and empty output Dataframe for complete CPD data.
    pp = CPD and TFBS is on + strand
    pm = CPD is on + strand and TFBS is on - strand
    mp = CPD is on - strand and TFBS is on + strand
    mm = CPD and TFBS is on - strand
    """
    all_same_seqs = pd.DataFrame()
    all_opp_seqs = pd.DataFrame()

    print('Counting CPDs...')
    for c in progressbar.progressbar(chroms):
        p_damages_chr = p_damages[p_damages["chr"] == c]
        m_damages_chr = m_damages[m_damages["chr"] == c]
        p_windows = sites.loc[(sites["chr"] == c) & (sites['strand'] == '+'), ['start', 'end', 'seq']]
        m_windows = sites.loc[(sites["chr"] == c) & (sites['strand'] == '-'), ['start', 'end', 'seq']]

        p_same_seqs = count_it(c, p_windows, window_sz,
                                    m_damages_chr[["start", "counts", "seq"]].to_records(index=False), True)
        m_same_seqs = count_it(c, m_windows, window_sz,
                                    p_damages_chr[["start", "counts", "seq"]].to_records(index=False), False)  #False
        p_opp_seqs = count_it(c, p_windows, window_sz,
                                   p_damages_chr[["start", "counts", "seq"]].to_records(index=False), True)
        m_opp_seqs = count_it(c, m_windows, window_sz,
                                   m_damages_chr[["start", "counts", "seq"]].to_records(index=False), False)  #False

        def concat_seqs(s, a_df):
            s_df = pd.DataFrame(s)
            tot_df = pd.concat([a_df, s_df])
            return tot_df

        all_same_seqs = concat_seqs(p_same_seqs, all_same_seqs)
        all_same_seqs = concat_seqs(m_same_seqs, all_same_seqs)
        all_opp_seqs = concat_seqs(p_opp_seqs, all_opp_seqs)
        all_opp_seqs = concat_seqs(m_opp_seqs, all_opp_seqs)

    all_same_seqs.to_csv(f"{output_path}/{tf_file}_{exp}_{kmer}mer_dipy_raw_counts_same.csv")
    all_opp_seqs.to_csv(f"{output_path}/{tf_file}_{exp}_{kmer}mer_dipy_raw_counts_opp.csv")



if __name__ == "__main__":
    main_dir = sys.argv[1]
    genome_path = sys.argv[2]
    data_path = sys.argv[3]
    kmer = int(sys.argv[4])
    tf_window_size = int(sys.argv[5])
    experiments_str = sys.argv[6]
    tf = sys.argv[7]
    experiments = experiments_str.split(",")
    tf_path = f'{main_dir}/results/TFBS'
    output_path = f'{main_dir}/results/TF_damage_data'

    for exp in experiments:
        print(f"Processing {tf} experiment {exp}...")
        proc_TFBS_dams(f'{tf}_top_{tf_window_size}_seq', kmer, exp, data_path, output_path, tf_path, genome_path)
        count_damages(f'{tf}_top_{tf_window_size}_seq', kmer, exp, output_path, tf_path)