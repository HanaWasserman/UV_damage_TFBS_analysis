import sys
import os
import pandas as pd
pd.options.mode.chained_assignment = None
import progressbar
from utils.pipeline_functions import rev_complement, adjust_kmer, filter_kmer, init4mers, init6mers

def prep_background_regions(dhs_file, data_path, output_path, genome_path):
    '''
    Subset intergenic, open chromatin regions for background modeling
    '''

    gen_fa_fai, gen_fa = f'{genome_path}/hg19/hg19.fa.fai', f'{genome_path}/hg19/hg19.fa'

    os.system(f"bedtools getfasta -fi {gen_fa} -bed {data_path}/{dhs_file}.bed -fo {output_path}/{dhs_file}_tmp_seq -bedOut")
    dhs_bed = pd.read_table(f'{data_path}/{dhs_file}.bed', header=None, names=['chrom','start','end'])
    seqs = pd.read_table(f"{output_path}/{dhs_file}_tmp_seq", header=None)
    dhs_bed['seq'] = seqs[0].str.upper()
    dhs_bed_complete = rev_complement(dhs_bed)
    os.remove(f"{output_path}/{dhs_file}_tmp_seq")
    dhs_bed_reordered = dhs_bed_complete[['chrom','start','end','seq','seq_rv']].copy()
    dhs_bed_reordered.to_csv(f"{output_path}/{dhs_file}_seq.bed", header=None, index=False, sep='\t')

def proc_background_dams(dhs_file, kmer, exp, output_path, data_path, genome_path, offset=1):

    fasta = f"{genome_path}/hg19/hg19.fa"
    f_path = f"{output_path}/{dhs_file}_{exp}_{kmer}mer"

    def prep_damage(strand_txt):
        print(f'Preparing: {f_path}')
        os.system(f"bedtools intersect -u -a {data_path}/{exp}_{strand_txt}_agg_filt.bed -b {data_path}/{dhs_file}.bed > {f_path}_{strand_txt}strand.bed")  # OUTPUT1, OUTPUT2
        df = pd.read_table(f'{f_path}_{strand_txt}strand.bed', header=None, names=['chr', 'st', 'end', 'counts'])
        return df

    prep_df_plus = prep_damage('plus')
    prep_df_minus = prep_damage('minus')

    adj_df_plus = adjust_kmer(prep_df_plus, -(offset + 1), (offset + 2), f_path, 'plus', ' -s', fasta)
    adj_df_minus = adjust_kmer(prep_df_minus, -(offset + 1), (offset + 2), f_path, 'minus', '', fasta)

    print(f'Filtering: {exp}')
    filter_kmer(adj_df_plus, offset, 'plus', f'{output_path}/{dhs_file}_intersect_{exp}_{kmer}mer')
    filter_kmer(adj_df_minus, offset, 'minus', f'{output_path}/{dhs_file}_intersect_{exp}_{kmer}mer')

    f_path = f"{output_path}/{dhs_file}_{exp}_{kmer}mer"
    os.remove(f"{f_path}_minusstrand.bed")
    os.remove(f"{f_path}_plusstrand.bed")
    os.remove(f"{f_path}_minus_tmp_tmp_seqs")
    os.remove(f"{f_path}_plus_tmp_tmp_seqs")
    os.remove(f"{f_path}_tmp_tmp.bed")

def prep_background_dams(dhs_file, exp, kmer, output_path, offset=1):
    damage_plus = pd.read_table(f'{output_path}/{dhs_file}_intersect_{exp}_{kmer}mer_dipy_proc_plus.bed', names=['chr', 'start', 'end','damage_count','seq', 'strand'])
    damage_minus = pd.read_table(f'{output_path}/{dhs_file}_intersect_{exp}_{kmer}mer_dipy_proc_minus.bed', names=['chr', 'start', 'end', 'damage_count', 'seq', 'strand'])

    bk_bed = pd.read_table(f"{output_path}/{dhs_file}_seq.bed", header=None, names=['chr', 'start', 'end', 'seq', 'seq_rc'])
    bk_bed['seq'] = bk_bed['seq'].str.upper()
    bk_bed['seq_rc'] = bk_bed['seq_rc'].str.upper()
    bk_bed_f = bk_bed[['chr', 'start', 'seq']].values.tolist()
    bk_bed_r = bk_bed[['chr', 'end', 'seq_rc']].values.tolist()

    bk_bed_f_list = list()
    for site in progressbar.progressbar(bk_bed_f):
        c, st, sq = site
        seq_l = len(sq)
        for i in range(seq_l - (kmer-1)):
            seq_n = sq[i:i + kmer]
            if 'N' in seq_n:
                continue
            if seq_n[(offset+1):(offset+3)] in ['CT', 'CC', 'TC', 'TT']:
                bk_bed_f_list.append((c, st + i, st + i + kmer, seq_n))
    bk_bed_f_df = pd.DataFrame()
    bk_bed_f_df['chr'] = [ch[0] for ch in bk_bed_f_list]
    bk_bed_f_df['start'] = [st[1] for st in bk_bed_f_list]
    bk_bed_f_df['end'] = [en[2] for en in bk_bed_f_list]
    bk_bed_f_df['seq'] = [se[3] for se in bk_bed_f_list]
    bk_bed_f_df = bk_bed_f_df.drop_duplicates(subset=['chr', 'start'])

    bk_bed_r_list = list()
    for site in progressbar.progressbar(bk_bed_r):
        c, en, sq = site
        seq_l = len(sq)
        for i in range(seq_l - (kmer-1)):
            seq_n = sq[i:i + kmer]
            if 'N' in seq_n:
                continue
            if seq_n[(offset+1):(offset+3)] in ['CT', 'CC', 'TC', 'TT']:
                bk_bed_r_list.append((c, en - i - kmer, en - i, seq_n))

    bk_bed_r_df = pd.DataFrame()
    bk_bed_r_df['chr'] = [ch[0] for ch in bk_bed_r_list]
    bk_bed_r_df['start'] = [st[1] for st in bk_bed_r_list]
    bk_bed_r_df['end'] = [en[2] for en in bk_bed_r_list]
    bk_bed_r_df['seq'] = [se[3] for se in bk_bed_r_list]
    bk_bed_r_df = bk_bed_r_df.drop_duplicates(subset=['chr', 'start'])

    bk_bed_dam_f = pd.merge(bk_bed_f_df[['chr', 'start', 'seq']], damage_minus, on=['chr', 'start', 'seq'], how='outer').fillna(0)
    bk_bed_dam_f = bk_bed_dam_f.drop(['chr', 'start'], axis=1)
    bk_bed_r_df = bk_bed_r_df.rename(columns={'seq_rc': 'seq'})
    bk_bed_dam_r = pd.merge(bk_bed_r_df[['chr', 'start', 'seq']], damage_plus, on=['chr', 'start', 'seq'], how='outer').fillna(0)
    bk_bed_dam_r = bk_bed_dam_r.drop(['chr', 'start'], axis=1)
    bk_bed_dam = pd.concat([bk_bed_dam_f, bk_bed_dam_r])
    bk_bed_dam.to_pickle(f'{output_path}/{dhs_file}_intersect_{exp}_{kmer}mer_background_all_data.pkl')

def make_damage_model(dhs_file, kmer, exp, background_path):
    def get_stats(damage_dat):
        if kmer == 6:
            damage_kmers = init6mers()
        elif kmer == 4:
            damage_kmers = init4mers()
            damage_dat['new_seq'] = damage_dat['seq'].str[1:5]
            damage_dat = damage_dat.drop(columns=['seq'])
            damage_dat = damage_dat.rename(columns={'new_seq': 'seq'})

        grouped_data = damage_dat.groupby('seq')['damage_count']
        precomputed_stats = {
            f: {
                'sum': group.sum(),
                'count': group.count(),
                'mean': group.mean(),
                'var': group.var()
            }
            for f, group in grouped_data
        }
        stats_l = []
        for f in progressbar.progressbar(damage_kmers):
            if f in precomputed_stats:
                stats = precomputed_stats[f]
                kmer_n = stats['count']
                damage_total = stats['sum']
                m = stats['mean']
                v = stats['var']
                stats_l.append((f, kmer_n, damage_total, m, v))

        stats_df = pd.DataFrame(stats_l, columns=['seq', 'kmer_count', 'damage_count', 'mean', 'var'])
        return stats_df

    dat = pd.read_pickle(f'{background_path}/{dhs_file}_intersect_{exp}_{kmer}mer_background_all_data.pkl')
    stats_all = get_stats(dat)
    stats_all.to_csv(f'{background_path}/{dhs_file}_intersect_{exp}_initial_damageability_{kmer}mer_background.csv')


if __name__ == "__main__":
    main_dir = sys.argv[1]
    dhs_file = sys.argv[2]
    genome_path = sys.argv[3]
    kmer = int(sys.argv[4])
    data_path = sys.argv[5]
    experiments_str = sys.argv[6]
    experiments = experiments_str.split(",")
    output_path = f'{main_dir}/results/background'

    prep_background_regions(dhs_file, data_path, output_path, genome_path)
    for exp in experiments:
        proc_background_dams(dhs_file, kmer, exp, output_path, data_path, genome_path)
        prep_background_dams(dhs_file, exp, kmer, output_path)
        make_damage_model(dhs_file, kmer, exp, output_path)