import pandas as pd
import os
import pickle
from utils.pipeline_functions import bootstrap_repair_tetramer, bootstrap_repair_cumulative
import sys


def sim_repair(samp_n, tf, dhs_file, p, exp1, exp2, strand, tf_window_size, repair_window, background_path, tf_data_path, output_path, mode="cumulative"):
    print(f'position: {p}')

    with open(f'{background_path}/{dhs_file}_{exp1}_vs_{exp2}_kmer_repair_dict.pkl', 'rb') as fp:
        repair_dict = pickle.load(fp)
    def get_data(exp, cols):
        data_exp = pd.read_csv(f'{tf_data_path}/{tf}_top_{tf_window_size}_seq_{exp}_6mer_dipy_raw_counts_{strand}.csv',
                               names=['chr', 'start', 'pos', 'seq', 'count'], index_col=0, skiprows=1)
        data_exp['4mer'] = data_exp['seq'].str[1:5]
        data_exp['pos'] = data_exp['pos'] - (tf_window_size-1) + 2
        return data_exp[cols]

    cols = ['chr', 'start', 'pos', '4mer', 'count']
    data_exp_1 = get_data(exp1, cols)
    data_exp_2 = get_data(exp2, cols)
    data_exp_1_2 = pd.merge(data_exp_1[cols], data_exp_2[cols], on=['chr', 'start', 'pos', '4mer'],
                            suffixes=['_1', '_2'],
                            how='outer').fillna(0)

    if mode == "tetramer":
        pos_repair = bootstrap_repair_tetramer(p, data_exp_1_2, repair_dict, samp_n)
        pos_repair_df = pd.DataFrame(pos_repair, columns=['pos', '4mer', 'count_1', 'count_2', 'pval'])
    elif mode == 'cumulative':
        raw_repair_sim = bootstrap_repair_cumulative(p, data_exp_1_2, repair_dict, samp_n)

    #pos_repair_df.to_csv(f"{output_path}/{tf}_top_{dhs_file}_{exp1}_vs_{exp2}_sim_{samp_n}_pos{p}_{mode}_{strand}.csv", index=False, header=False)
    if raw_repair_sim is not None:
        raw_repair_sim['strand'] = strand
        raw_repair_sim['TF'] = tf
        raw_repair_sim.to_csv(f"{output_path}/raw_simulations/{tf}_top_{dhs_file}_{exp1}_vs_{exp2}_sim_{samp_n}_pos{p+repair_window}_{mode}_raw_{strand}.csv", index=False)
    else:
        return


if __name__ == "__main__":
    main_dir = sys.argv[1]
    tf = sys.argv[2]
    exp1 = sys.argv[3]
    exp2 = sys.argv[4]
    mode = sys.argv[5]
    strand = sys.argv[6]
    dhs_file = sys.argv[7]
    tf_window_size = sys.argv[8]
    samp_n = sys.argv[9]
    repair_window = sys.argv[10]
    output_path = f"{main_dir}/results/analysis/repair"
    background_path = f"{main_dir}/results/background"
    tf_data_path = f"{main_dir}/results/tf_damage_data"
    pos = int(os.getenv("SLURM_ARRAY_TASK_ID"))
    p = pos - repair_window


    sim_repair(samp_n, tf, dhs_file, pos, exp1, exp2, strand, tf_window_size, repair_window, background_path, tf_data_path, output_path)
