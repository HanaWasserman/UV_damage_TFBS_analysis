import pandas as pd
import sys
import pickle
import functools as ft
from utils.pipeline_functions import init4mers

def make_repair_dict(dhs_file, background_path, exp1, exp2):
    exp_1_same = pd.read_table(f'{background_path}/{dhs_file}_intersect_{exp1}_6mer_dipy_proc_plus.bed',
                               usecols=[0, 1, 3, 4], names=['chr', 'start', 'count', 'seq'])
    exp_1_opp = pd.read_table(f'{background_path}/{dhs_file}_intersect_{exp1}_6mer_dipy_proc_minus.bed',
                              usecols=[0, 1, 3, 4], names=['chr', 'start', 'count', 'seq'])
    exp_2_same = pd.read_table(f'{background_path}/{dhs_file}_intersect_{exp2}_6mer_dipy_proc_plus.bed',
                               usecols=[0, 1, 3, 4], names=['chr', 'start', 'count', 'seq'])
    exp_2_opp = pd.read_table(f'{background_path}/{dhs_file}_intersect_{exp2}_6mer_dipy_proc_minus.bed',
                              usecols=[0, 1, 3, 4], names=['chr', 'start', 'count', 'seq'])
    exp_1_same['4mer'] = exp_1_same['seq'].str[1:5]
    exp_1_opp['4mer'] = exp_1_opp['seq'].str[1:5]
    exp_2_same['4mer'] = exp_2_same['seq'].str[1:5]
    exp_2_opp['4mer'] = exp_2_opp['seq'].str[1:5]

    cols = ['chr', 'start', 'count', '4mer']
    exp_same = [exp_1_same[cols], exp_2_same[cols]]
    exp_same_merge = ft.reduce(
        lambda left, right: pd.merge(left, right, on=['chr', 'start', '4mer'], how='outer', suffixes=['_1', '_2']),
        exp_same).fillna(0)
    exp_opp = [exp_1_opp[cols], exp_2_opp[cols]]
    exp_opp_merge = ft.reduce(
        lambda left, right: pd.merge(left, right, on=['chr', 'start', '4mer'], how='outer', suffixes=['_1', '_2']),
        exp_opp).fillna(0)
    exp = pd.concat([exp_same_merge, exp_opp_merge])

    repair_kmer_dict = init4mers()
    for kmer in repair_kmer_dict:
        repair_kmer_dict[kmer] = exp.loc[exp['4mer'] == kmer, ['count_1', 'count_2']]

    print(f'Saving dictionary for {exp1} vs. {exp2}...')
    with open(f'{background_path}/{dhs_file}_{exp1}_vs_{exp2}_kmer_repair_dict.pkl', 'wb') as fp:
        pickle.dump(repair_kmer_dict, fp)


if __name__ == "__main__":
    main_dir = sys.argv[1]
    dhs_file = sys.argv[2]
    kmer = int(sys.argv[3])
    exp1 = sys.argv[4]
    exp2 = sys.argv[5]
    background_path = f'{main_dir}/results/background'

    make_repair_dict(dhs_file, background_path, exp1, exp2)
