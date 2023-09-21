import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(font_scale=1.5)
sns.set_style("whitegrid")

path = "/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/snakemake_assembly/binning/metawrap/BIN_REFINEMENT"

sample_name_mapping = {
    "BGG_BOSW": "0_NA_B",
    "BGG_POSDC": "1_11_P",
    "BGG_KOSDC": "1_16_K",
    "BGG_LOSDA": "2_12_L",
    "BGG_NOSDA": "2_16_N",
    "BGG_ROSDA": "2_19_R",
    "BGG_QOSDA": "3_06_Q",
    "BGG_MOSDA": "3_16_M",
    "BGG_OOSDC": "3_19_O",
}

# Order based on the sample_name_mapping
order = list(sample_name_mapping.values())

# Dictionary to store data from each file
data_dict = {}

for folder in os.listdir(path):
    if folder.startswith("BGG"):
        for file in os.listdir(os.path.join(path, folder)):
            if file == "metawrap_50_10_bins.stats":
                sample_name = folder.split("-")[0]
                sample_name = sample_name_mapping[sample_name]
                print(f"----------File: {sample_name}----------")
                data_dict[sample_name] = pd.read_csv(os.path.join(path, folder, file), sep="\t")

data_concat = pd.concat(data_dict.values(), keys=data_dict.keys()).reset_index()
data_concat = data_concat.rename(columns={"level_0": "sample"})

print(data_concat)

# Plot for completeness
plt.figure(figsize=(14, 6))
sns.catplot(
    data=data_concat, 
    x='sample', 
    y='completeness', 
    kind='box', 
    height=5, 
    aspect=1.5, 
    dodge=True,
    showfliers=True,
    order=order
)

sns.stripplot(
    data=data_concat, 
    x='sample', 
    y='completeness', 
    jitter=True, 
    color='lightgrey', 
    dodge=True, 
    alpha=0.3
)

plt.xticks(rotation=45)
plt.tight_layout()
plt.xlabel("Sample")
plt.ylabel("Completeness [%]")
plt.savefig("../output/MAG_completeness_boxplot.png", dpi=300)

# Plot for contamination
plt.figure(figsize=(14, 6))
sns.catplot(
    data=data_concat, 
    x='sample', 
    y='contamination', 
    kind='box', 
    height=5, 
    aspect=1.5, 
    dodge=True,
    showfliers=True,
    order=order
)

sns.stripplot(
    data=data_concat, 
    x='sample', 
    y='contamination', 
    jitter=True, 
    color='lightgrey', 
    dodge=True, 
    alpha=0.3
)

plt.xticks(rotation=45)
plt.tight_layout()
plt.xlabel("Sample")
plt.ylabel("Contamination [%]")
plt.savefig("../output/MAG_contamination_boxplot.png", dpi=300)