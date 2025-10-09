import pandas as pd
import os
import glob
import re
import sys


def analyze_repair_sims(tf, tf_len, exp1, exp2, samp_n, dhs_file, analysis_path, mode="cumulative"):

    def analyze_strand(strand):
        def get_scaling_factor(df, start, stop):
            flank = df[df['pos'].between(start, stop)]
            return (flank['count_2'] / flank['count_2_sim']).mean()

        csv_files = glob.glob(
            f"{analysis_path}/raw_simulations/{tf}_top_{dhs_file}_{exp1}_vs_{exp2}_sim_{samp_n}_pos*_{mode}_raw_{strand}.csv")
        dfs = []
        print(f"{analysis_path}/raw_simulations/{tf}_top_{dhs_file}_{exp1}_vs_{exp2}_sim_{samp_n}_pos*_{mode}_raw_{strand}.csv")
        for file in csv_files:
            match = re.search(r'pos(\d+)', os.path.basename(file))
            if match:
                pos_number = int(match.group(1))
            else:
                pos_number = None
            df = pd.read_csv(file)
            df['pos'] = pos_number
            dfs.append(df)

        sim = pd.concat(dfs, ignore_index=True)
        sim['pos'] = sim['pos'] - 30
        print(sim)
        sim_agg = sim.groupby(by=['pos', 'count_1', 'count_2', 'strand','TF'], as_index=False).mean()
        sim_agg = sim_agg.sort_values(by='pos')
        scaling_l = get_scaling_factor(sim_agg, sim_agg['pos'].min(), -((tf_len//2)+10))
        scaling_r = get_scaling_factor(sim_agg, (tf_len // 2) + 10, sim_agg['pos'].max())
        scaling_0 = (scaling_l+scaling_r)/2


        sim['sim_count_2_scaled'] = float(0)
        sim['scaling_factor'] = float(0)
        sim.loc[sim['pos'] < 0, 'count_2_sim_scaled'] = sim.loc[sim['pos'] < 0, 'count_2_sim'] * scaling_l
        sim.loc[sim['pos'] > 0, 'count_2_sim_scaled'] = sim.loc[sim['pos'] > 0, 'count_2_sim'] * scaling_r
        sim.loc[sim['pos'] == 0, 'count_2_sim_scaled'] = sim.loc[sim['pos'] == 0, 'count_2_sim'] * scaling_0
        sim.loc[sim['pos'] < 0, 'scaling_factor'] = scaling_l
        sim.loc[sim['pos'] > 0, 'scaling_factor'] = scaling_r
        sim.loc[sim['pos'] == 0, 'scaling_factor'] = scaling_0

        repair_res = []
        for p in list(sim['pos'].unique()):
            pos_sim = sim.loc[sim['pos'] == p]
            pos_count_1 = pos_sim['count_1'].max()
            pos_count_2 = pos_sim['count_2'].max()

            pval_hi = pos_sim[pos_sim['count_2_sim'] > pos_count_2].shape[0] / pos_sim.shape[0]
            pval_hi_scaled = pos_sim[pos_sim['count_2_sim_scaled'] > pos_count_2].shape[0] / pos_sim.shape[0]
            pval_lo = pos_sim[pos_sim['count_2_sim'] < pos_count_2].shape[0] / pos_sim.shape[0]
            pval_lo_scaled = pos_sim[pos_sim['count_2_sim_scaled'] < pos_count_2].shape[0] / pos_sim.shape[0]

            mean = pos_sim['count_2_sim'].mean()
            mean_scaled = pos_sim['count_2_sim_scaled'].mean()
            var = pos_sim['count_2_sim'].var()
            var_scaled = pos_sim['count_2_sim_scaled'].var()
            std = pos_sim['count_2_sim'].std()
            std_scaled = pos_sim['count_2_sim_scaled'].std()

            repair_res.append(
                (p, pos_count_1, pos_count_2, mean_scaled, var_scaled, std_scaled, pval_hi_scaled, pval_lo_scaled, mean, var, std, pval_hi, pval_lo))

        repair_res_df = pd.DataFrame(repair_res,
                                  columns=['pos', 'count_1', 'count_2', 'mean_sc', 'var_sc', 'std_sc', 'pval_hi_sc', 'pval_lo_sc',
                                           'mean', 'var', 'std', 'pval_hi', 'pval_lo'])
        repair_res_df["TF"] = tf
        repair_res_df["strand"] = strand
        repair_res_df["scal_factor_l"] = scaling_l
        repair_res_df["scal_factor_r"] = scaling_r
        return repair_res_df

    repair_same = analyze_strand("same")
    repair_opp = analyze_strand("opp")
    repair_all = pd.concat([repair_same, repair_opp])
    repair_all.to_csv(f"{analysis_path}/{tf}_top_{dhs_file}_{exp1}_vs_{exp2}_{samp_n}_sim_scaled_{mode}.csv",
                  index=False, header=False)


if __name__ == "__main__":
    main_dir = sys.argv[1]
    exp1 = sys.argv[2]
    exp2 = sys.argv[3]
    dhs_file = sys.argv[4]
    samp_n = sys.argv[5]
    tf_len = int(sys.argv[6])
    tf = sys.argv[7]
    analysis_path = f"{main_dir}/results/analysis/repair"


    analyze_repair_sims(tf, tf_len, exp1, exp2, samp_n, dhs_file, analysis_path)

