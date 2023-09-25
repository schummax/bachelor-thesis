import pandas as pd
import os

path = "/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/snakemake_assembly/stats/bin_quality"

pass_list = []
rows_list = []
for file in os.listdir(path):
    pass_count = 0
    df = pd.read_csv(os.path.join(path, file), sep="\t")
    for index, row in df.iterrows():
        if row['pass.GUNC'] == True:
            pass_count += 1
    pass_list.append(pass_count)
    rows_list.append(df.shape[0])
    print(f"{file} -> {pass_count}/{df.shape[0]}")

total_pass = sum(pass_list)
total_rows = sum(rows_list)
print(f"Total pass: {total_pass}/{total_rows}")
print(f"Percentage: {round(total_pass/total_rows*100, 2)}%")