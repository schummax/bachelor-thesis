import pandas as pd
import os

gtdbtk_bac = pd.read_csv(
    "/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/snakemake_assembly/stats/gtdbtk/gtdbtk.Rigou2022_bac120.summary.tsv",
    sep="\t",
)
gtdbtk_ar = pd.read_csv(
    "/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/snakemake_assembly/stats/gtdbtk/gtdbtk.Rigou2022_ar53.summary.tsv",
    sep="\t",
)

gtdbtk = pd.concat([gtdbtk_bac, gtdbtk_ar])

NODE_143 = "Polyangiales"
NODE_402 = "Chlamydiales"
M = "Mycolicibacterium"

NODE_143_count = 0
NODE_402_count = 0
M_count = 0

for index, row in gtdbtk.iterrows():
    if NODE_143 in row["classification"] or NODE_402 in row["classification"] or M in row["classification"]:
        if NODE_143 in row["classification"]:
            print("Found Polyangiales") 
            NODE_143_count += 1
        elif NODE_402 in row["classification"]:
            print("Found Chlamydiales")
            NODE_402_count += 1
        elif M in row["classification"]:
            print("Found Mycolicibacterium")
            M_count += 1
        else:
            print("Error")
            
print(f"Polyangiales: {NODE_143_count}")
print(f"Chlamydiales: {NODE_402_count}")
print(f"Mycolicibacterium: {M_count}")