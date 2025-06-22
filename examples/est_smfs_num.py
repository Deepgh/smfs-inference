

if __name__ == "__main__":
    import os
    import pandas as pd
    from pathlib import Path
    from smfs.estimation.smfs_fit_num import compute
    path = '/project/yuvalsim/Deep/project2/human_data/sim_hum_data/final_approach/smp_mu_diff_site/seeded'
    data_file = 'edited_sim_data_smp_mu_diff_site.csv'
    unq_mut_sites_path = path
    unq_mut_sites_file = 'unq_mut_sites_mu_1000_sim_data_diff_site.csv'
    ac_limit = 1000
    output_path = path
    output_filename = f'smfs_num_{ac_limit}_sim_data_mu_diff_site.csv'
    Path(output_path).mkdir(parents=True, exist_ok=True)
    data_df = pd.read_csv(os.path.join(path, data_file))
    pmf_data_file = os.path.join(unq_mut_sites_path, unq_mut_sites_file)
    df_pmf = pd.read_csv(pmf_data_file)
    f = compute(data_df['AC'].to_numpy(), df_pmf['Unq muts sites'].to_numpy(), ac_limit)
    new_df = pd.DataFrame()
    for i in range(len(f)):
        data = pd.DataFrame([{'Count': i+1, 'pred prob':f[i]}])
        new_df = pd.concat([new_df, data], ignore_index=True)
    new_df.to_csv(os.path.join(output_path, output_filename), index=False)