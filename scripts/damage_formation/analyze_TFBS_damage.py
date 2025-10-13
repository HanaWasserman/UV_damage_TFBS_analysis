import sys
import os
import pandas as pd
import numpy as np
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
from utils.pipeline_functions import fill_pos_gaps

def predict_TFBS_damage(tf_file, exp, kmer, tf_window_size, dhs_file, background_path, output_path, tf_path):

    kmer_bk = pd.read_csv(f'{background_path}/{dhs_file}_intersect_{exp}_initial_damageability_{kmer}mer_background.csv', index_col=0)

    def pred_damage_poisson(cnts):
        damage_pred = pd.merge(cnts, kmer_bk[['seq', 'mean', 'var']], on='seq')
        damage_pred['pred_mean'] = damage_pred['count'] * damage_pred['mean']
        damage_pred['pred_var'] = damage_pred['count'] * damage_pred['var']
        return damage_pred

    def aggregate_damage(kmer_counts, strand):
        kmer_counts_agg = kmer_counts.groupby(by='pos', as_index=False).sum().reset_index(drop=True)
        kmer_counts_agg['std_dev'] = np.sqrt(kmer_counts_agg['pred_mean'])

        kmer_counts_agg = fill_pos_gaps(kmer_counts_agg)
        kmer_counts_agg[['pos', 'count', 'mean', 'var', 'pred_mean', 'pred_var', 'std_dev']].to_csv(
            f"{output_path}/{tf_file}_{exp}_{kmer}mer_predicted_damage_{dhs_file}_ROI_{strand}.csv", index=False)

    print(f'Predicting damages for {tf_file}')
    kmers_plus = pd.read_csv(f"{tf_path}/{tf_file}_{kmer}mers_ROI_plus.csv")
    kmers_minus = pd.read_csv(f"{tf_path}/{tf_file}_{kmer}mers_ROI_minus.csv")
    kmers_plus['pos'] = kmers_plus['pos'] - (tf_window_size-1)
    kmers_minus['pos'] = kmers_minus['pos'] - (tf_window_size-1)

    p_kmer_counts = pred_damage_poisson(kmers_plus.loc[kmers_plus['pos'].between(-(tf_window_size/2), tf_window_size/2)])
    m_kmer_counts = pred_damage_poisson(kmers_minus.loc[kmers_minus['pos'].between(-(tf_window_size/2), tf_window_size/2)])
    p_kmer_counts.to_csv(f"{output_path}/{tf_file}_{exp}_{kmer}mer_predicted_damage_{dhs_file}_ROI_plus_ALL.csv", index=False)
    m_kmer_counts.to_csv(f"{output_path}/{tf_file}_{exp}_{kmer}mer_predicted_damage_{dhs_file}_ROI_minus_ALL.csv", index=False)
    aggregate_damage(p_kmer_counts, 'same')
    aggregate_damage(m_kmer_counts, 'opp')

def analyze_TFBS_damage(tf, tf_file, exp, kmer, tf_window_size, anal_window, tfbs_buffer, dhs_file, tf_data_path, output_path):

    if not os.path.exists(output_path):
        os.system(f"mkdir {output_path}")
    def get_data(strand):
        d = pd.read_csv(f'{tf_data_path}/{tf_file}_seq_{exp}_6mer_dipy_raw_counts_{strand}.csv',
            index_col=0, names=['chr', 'start', 'pos','seq','count'], skiprows=1)
        d_agg = d[['pos', 'count']].groupby(by='pos', as_index=False).sum().reset_index(drop=True)
        d_agg = fill_pos_gaps(d_agg)
        if kmer==6:
            d_agg['pos'] = d_agg['pos'] - (tf_window_size-1) + 2.5
        d_agg = d_agg.loc[d_agg['pos'].between(-anal_window, anal_window)].copy().reset_index(drop=True)
        return d_agg

    d_opp = get_data('opp')
    d_same = get_data('same')

    def get_pred(strand):
        p = pd.read_csv(f'{tf_data_path}/{tf_file}_{exp}_{kmer}mer_predicted_damage_{dhs_file}_ROI_{strand}.csv')
        p_subset = p.loc[p['pos'].between(-anal_window, anal_window)].copy()
        p['pos'] = p['pos'] + 2.5
        p_subset['pos'] = p_subset['pos'] + 2.5
        return p, p_subset

    p_opp_all, p_opp = get_pred('opp')
    p_same_all, p_same = get_pred('same')

    def calc_pval(df_row):
        if df_row['count'] > df_row['pred_mean_scaled']:
            pval = 1 - poisson.cdf(df_row['count'], df_row['pred_mean_scaled'])
            enriched = True
        else:
            pval = poisson.cdf(df_row['count'], df_row['pred_mean_scaled'])
            enriched = False
        return pd.Series([pval, enriched], index=['pval', 'enriched'])

    def calc_zscore(df_row):
        if df_row['std_dev'] != 0:
            zscore = (df_row['count'] - df_row['pred_mean']) / df_row['std_dev']
            zscore_scal = (df_row['count'] - df_row['pred_mean_scaled']) / df_row['std_dev_scal']
            return pd.Series([zscore, zscore_scal], index=['zscore', 'zscore_scal'])
        else:
            return pd.Series([0,0], index=['zscore', 'zscore_scal'])

    def merge_dat_pred(d, p):
        d_p_merge = pd.merge(d, p[['pos', 'pred_mean', 'pred_var', 'std_dev']], on='pos')
        scaling_factor_l = d_p_merge.loc[d_p_merge['pos'].between(-tfbs_buffer-15,-tfbs_buffer), 'count'].mean() / d_p_merge.loc[d_p_merge['pos'].between(-tfbs_buffer-15,-tfbs_buffer), 'pred_mean'].mean()
        scaling_factor_r = d_p_merge.loc[d_p_merge['pos'].between(tfbs_buffer, tfbs_buffer+15), 'count'].mean() / d_p_merge.loc[d_p_merge['pos'].between(tfbs_buffer, tfbs_buffer+15), 'pred_mean'].mean()
        d_p_merge['pred_mean_scaled'] = float(0)
        d_p_merge.loc[d_p_merge['pos'].between(-anal_window, -.5),'pred_mean_scaled'] = d_p_merge.loc[d_p_merge['pos'].between(-anal_window, -.5),'pred_mean']*scaling_factor_l
        d_p_merge.loc[d_p_merge['pos'].between(1.5, anal_window), 'pred_mean_scaled'] = d_p_merge.loc[d_p_merge['pos'].between(1.5,anal_window), 'pred_mean'] * scaling_factor_r
        d_p_merge.loc[d_p_merge['pos'].between(-anal_window, -.5), 'pred_var_scaled'] = d_p_merge.loc[d_p_merge['pos'].between(-anal_window, -.5), 'pred_var'] * scaling_factor_l
        d_p_merge.loc[d_p_merge['pos'].between(1.5, anal_window), 'pred_var_scaled'] = d_p_merge.loc[d_p_merge['pos'].between(1.5, anal_window), 'pred_var'] * scaling_factor_r
        d_p_merge.loc[d_p_merge['pos']==0.5, 'pred_mean_scaled'] = d_p_merge.loc[d_p_merge['pos']==0.5, 'pred_mean'] * ((scaling_factor_r+scaling_factor_l)/2)
        d_p_merge.loc[d_p_merge['pos'] == 0.5, 'pred_var_scaled'] = d_p_merge.loc[d_p_merge['pos'] == 0.5, 'pred_var'] * ((scaling_factor_r + scaling_factor_l) / 2)
        d_p_merge['std_dev_scal'] = np.sqrt(d_p_merge['pred_mean_scaled'])
        d_p_merge[['pval', 'enriched']] = d_p_merge.apply(calc_pval, axis=1)
        d_p_merge[['zscore','zscore_scal']] = d_p_merge.apply(calc_zscore, axis=1)
        d_p_merge['scal_factor_l'] = scaling_factor_l
        d_p_merge['scal_factor_r'] = scaling_factor_r
        return d_p_merge.sort_values(by='pos'), scaling_factor_l, scaling_factor_r

    d_p_same, scal_factor_same_l, scal_factor_same_r = merge_dat_pred(d_same, p_same)
    d_p_opp, scal_factor_opp_l, scal_factor_opp_r = merge_dat_pred(d_opp, p_opp)
    d_p_same['strand'] = 'same'
    d_p_opp['strand'] = 'opp'

    d_p = pd.concat([d_p_same, d_p_opp])
    d_p_sorted = d_p.sort_values(by='pval')
    d_p_sorted_corrected_pvals = multipletests(d_p_sorted['pval'], method='fdr_bh')
    d_p_sorted['qval'] = d_p_sorted_corrected_pvals[1]

    cols = ['pos', 'count', 'pred_mean', 'std_dev', 'pred_mean_scaled', 'std_dev_scal', 'pval', 'qval', 'enriched', 'zscore', 'zscore_scal', 'scal_factor_l', 'scal_factor_r','strand']
    pv = d_p_sorted[cols]
    pv['TF'] = tf
    pv.to_csv(f'{output_path}/{tf_file}_{exp}_{dhs_file}_{kmer}mer_pvals_corrected.csv', header=False, index=False)

if __name__ == "__main__":
    main_dir = sys.argv[1]
    kmer = int(sys.argv[2])
    tf_window_size = int(sys.argv[3])
    analysis_window_size = int(sys.argv[4])
    dhs_file = sys.argv[5]
    experiments_str = sys.argv[6]
    tf_len = int(sys.argv[7])
    tf = sys.argv[8]
    experiments = experiments_str.split(",")
    background_path = f'{main_dir}/results/background'
    tf_data_path = f'{main_dir}/results/TF_damage_data'
    tf_path = f'{main_dir}/results/TFBS'
    analysis_path = f'{main_dir}/results/analysis/damage_formation'

    for exp in experiments:
        print(f"Analyzing {tf} experiment {exp}...")
        predict_TFBS_damage(f'{tf}_top_{tf_window_size}', exp, kmer, tf_window_size, dhs_file, background_path, tf_data_path, tf_path)
        analyze_TFBS_damage(tf, f'{tf}_top_{tf_window_size}', exp, kmer, tf_window_size, (tf_len//2+analysis_window_size), (tf_len//2+5), dhs_file, tf_data_path, analysis_path)