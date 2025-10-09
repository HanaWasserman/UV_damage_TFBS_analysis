import pandas as pd
import os
import numpy as np
import progressbar
from Bio.Seq import Seq
import logomaker as lm
from tqdm import tqdm

def init4mers():
    yys = ('CC', 'CT', 'TC', 'TT')
    ns = ('A', 'C', 'T', 'G')
    qn_l = list()
    for n1 in ns:
        for n2 in ns:
            for yy in yys:
                    qn_l.append((f'{n1}{yy}{n2}'))
    qn = dict.fromkeys(qn_l, 0)
    for q in qn:
        qn[q] = []
    return qn

def init6mers():
    yys = ('CC', 'CT', 'TC', 'TT')
    ns = ('A', 'C', 'T', 'G')
    qn_l = list()
    for n1 in ns:
        for n2 in ns:
            for n3 in ns:
                for n4 in ns:
                    for yy in yys:
                        qn_l.append((f'{n1}{n2}{yy}{n3}{n4}'))
    qn = dict.fromkeys(qn_l, 0)
    for q in qn:
        qn[q] = []
    return qn

def rev_complement(df):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    df['seq_rv'] = ""
    for i in progressbar.progressbar(range(df.shape[0])):
        df['seq_rv'][i] = "".join(complement.get(base, base) for base in reversed(df['seq'][i]))
    return df

def get_tf(bed_file, cols):
    bs = pd.read_table(f'{bed_file}.bed', header=None, usecols=[0, 1, 2, 3, 4, 5], names=['chr', 'start', 'end', 'seq', 'rv_seq', 'strand'])
    bs[["start", "end"]] = bs[["start", "end"]].to_numpy(dtype='int')
    ws = bs["end"][0] - bs["start"][0]
    return bs[cols], ws

def adjust_kmer(df, shift_left, shift_right, f_path, strand, strandedness, fasta):
    """
    Shifts CPD signal intervals to accurately map dinucleotide fasta sequence to each CPD using BedTools
    """
    df["strand"] = "-"
    df["name"] = ""
    df["st"] += shift_left
    df["end"] += shift_right
    df_ = df[['chr','st','end','counts','name','strand']].copy()
    df_.to_csv(f'{f_path}_tmp_tmp.bed', header=None, index=False, sep='\t')
    os.system(f"bedtools getfasta -fi {fasta} -bed {f_path}_tmp_tmp.bed -fo {f_path}_{strand}_tmp_tmp_seqs{strandedness} -bedOut")
    seqs = pd.read_table(f"{f_path}_{strand}_tmp_tmp_seqs", header=None)
    df_["seq"] = seqs[0].str.upper()
    if strand == 'plus':
        df_["ori"] = "+"
    return df_[['chr','st','end','counts','seq','strand']]

def filter_kmer(df, offset, strand, f_path):
    """
    Verifies CPD signal maps to dipyrimdine sequences and filters out any non-matches
    """
    all_n = df.shape[0]
    df["seq"] = df["seq"].str.upper()
    df_dipy = df[(df["seq"].str[(offset+1):(offset+3)] == "CC") | (df["seq"].str[(offset+1):(offset+3)] == "CT") | (df["seq"].str[(offset+1):(offset+3)] == "TC") | (df["seq"].str[(offset+1):(offset+3)] == "TT")]
    df_dipy_ = df_dipy[["chr", "st", "end", "counts", "seq", "strand"]].copy()
    dipy_n = str(df_dipy_.shape[0] / all_n)
    df_dipy_.to_csv(f"{f_path}_dipy_proc_{strand}.bed", header=None, index=False, sep='\t')
    print(f"{dipy_n}% dipyrimidines")
    return df_dipy_

def fill_pos_gaps(dat):
    new_pos = pd.DataFrame()
    new_pos['pos'] = np.arange(dat['pos'].min(), dat['pos'].max())
    return pd.merge(new_pos, dat, on=['pos'], how='outer').fillna(0)

def bootstrap_repair_cumulative(p, exp_dict, nyyn_dict, samp_n):
    """
    Bootstrap sampling of CPD repair dictionary to assess aggregate CPD counts at
    second timepoint, based on NYYN composition and starting damage levels
    :param p: position
    :param exp_dict: damages merged across both timepoints by genomic locus
    :param nyyn_dict: background tetranucleotide repair dictionary for sampling
    :param samp_n: bootstrap samples
    """
    exp_dict_pos = exp_dict.loc[exp_dict['pos'] == p]
    start_signal_pos = exp_dict_pos['count_1'].sum()
    if start_signal_pos == 0:
        #return pd.DataFrame([[p, 0, 0, 1, 1, 0, 0, 0]], columns=['pos','count_1','count_2','pval', 'pval_2', 'sim_mean','sim_var','sim_std']), None
        return None
    end_signal_pos = exp_dict_pos['count_2'].sum()
    sim_pos = list()
    for s in progressbar.progressbar(nyyn_dict):
        exp_dict_pos_seq = exp_dict_pos.loc[exp_dict_pos['4mer'] == s]
        start_signal_frm = exp_dict_pos_seq['count_1'].sum()
        if start_signal_frm == 0:
            continue
        nyyn_dict_seq = nyyn_dict[s].copy()
        while nyyn_dict_seq['count_1'].sum() < start_signal_frm:
            nyyn_dict_seq = pd.concat([nyyn_dict_seq, nyyn_dict_seq])
        nyyn_dict_seq_n = nyyn_dict_seq.shape[0]
        for i in range(samp_n):
            while True:
                samp = nyyn_dict_seq.sample(n=nyyn_dict_seq_n)
                sim_count_1_try = samp[samp['count_1'].cumsum() <= start_signal_frm]['count_1'].sum()
                if sim_count_1_try < start_signal_frm:
                    continue
                else:
                    samp_cum = samp[samp['count_1'].cumsum() <= start_signal_frm]
                    while True:
                        if samp_cum.shape[0]==1:
                            break
                        if samp_cum[:-1]['count_1'].sum() == start_signal_frm:
                            samp_cum = samp_cum[:-1].copy()
                            continue
                        else:
                            break
                    sim_count_2 = samp_cum['count_2'].sum()
                    sim_pos.append((s, sim_count_2))
                    break
    sim_pos_df = pd.DataFrame(sim_pos, columns=['4mer','count_2_frm_sim'])
    sim_pos_by_seq = {k: v.reset_index(drop=True) for k, v in sim_pos_df.groupby("4mer")}
    sim_pos_agg_df = pd.DataFrame({"count_2_sim": pd.concat([df["count_2_frm_sim"].reset_index(drop=True) for df in sim_pos_by_seq.values()], axis=1).sum(axis=1)})
    sim_pos_agg_df['pos'] = p
    sim_pos_agg_df['count_1'] = start_signal_pos
    sim_pos_agg_df['count_2'] = end_signal_pos
    #sim_pos_stats = []
    #pval_top = sim_pos_agg_df[sim_pos_agg_df['count_2_sim']>end_signal_pos].shape[0]/samp_n
    #pval_bottom = sim_pos_agg_df[sim_pos_agg_df['count_2_sim']<end_signal_pos].shape[0] / samp_n
    #sim_pos_stats.append((p, start_signal_pos, end_signal_pos, sim_pos_agg_df['count_2_sim'].mean(), sim_pos_agg_df['count_2_sim'].var(), sim_pos_agg_df['count_2_sim'].std(), pval_top, pval_bottom))
    #sim_pos_stats_df = pd.DataFrame(sim_pos_stats, columns=["pos","count_1","count_2","count_2_sim_mean","count_2_sim_var","count_2_sim_std","pval_top","pval_bottom"])
    return sim_pos_agg_df[["pos","count_1","count_2","count_2_sim"]]

def simulate_repair_for_tetramer(p, f_, obs1_obs2_pos, frm, samp_n):
    o_s = obs1_obs2_pos.loc[obs1_obs2_pos['4mer'] == f_]
    start_signal = o_s['count_1'].sum()
    if start_signal == 0:
        return None
    end_signal = o_s['count_2'].sum()
    shuff_o2_all = []
    for i in range(samp_n):
        np.random.seed(i)
        samp_shuff = frm[f_].sample(n=frm[f_].shape[0], replace=True)
        while samp_shuff['count_1'].sum() < start_signal:
            samp_shuff = pd.concat([samp_shuff, samp_shuff])
        shuff_o2 = samp_shuff[samp_shuff['count_1'].cumsum() <= start_signal]['count_2'].sum()
        shuff_o2_all.append(shuff_o2)
    pval = 1 - (sum(ele < end_signal for ele in shuff_o2_all) / len(shuff_o2_all))
    return (p, f_, start_signal, end_signal, pval)

def bootstrap_repair_tetramer(p, o1_o2, frm):
    results_by_pos = []
    obs1_obs2_pos = o1_o2.loc[o1_o2['pos'] == p]
    for f_ in tqdm(frm):
        result = simulate_repair_for_tetramer(p, f_, obs1_obs2_pos, frm)
        if result:
            results_by_pos.append(result)
    return results_by_pos