import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

gtdbtk_bac = pd.read_csv(
    "/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/snakemake_assembly/stats/gtdbtk/gtdbtk.Rigou2022_bac120.summary.tsv",
    sep="\t",
)
gtdbtk_ar = pd.read_csv(
    "/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/snakemake_assembly/stats/gtdbtk/gtdbtk.Rigou2022_ar53.summary.tsv",
    sep="\t",
)

gtdbtk = pd.concat([gtdbtk_bac, gtdbtk_ar])

# Extract taxonomic-level from the classification column
gtdbtk["genus"] = gtdbtk["classification"].str.split(";").str[5]
gtdbtk["class"] = gtdbtk["classification"].str.split(";").str[2]
gtdbtk["order"] = gtdbtk["classification"].str.split(";").str[3]

gtdbtk["genus"] = gtdbtk["genus"].str.replace("g__", "")
gtdbtk["class"] = gtdbtk["class"].str.replace("c__", "")
gtdbtk["order"] = gtdbtk["order"].str.replace("o__", "")

gtdbtk["order"] = gtdbtk["order"].replace("", "unclassified")

# Count the occurrence of each level and select the top 20
top_20_genera = gtdbtk["genus"].value_counts().head(20)
top_20_classes = gtdbtk["class"].value_counts().head(20)
top_20_orders = gtdbtk["order"].value_counts().head(20)

print(top_20_genera)
print(top_20_classes)
print(top_20_orders)

# Plot the top 20 genera
sns.set_theme(style="whitegrid")

plt.figure(figsize=(10, 8))
sns.barplot(x=top_20_genera.values, y=top_20_genera.index, color="#96B0CA")

plt.xlabel("Number of MAGs", fontsize=14)
plt.ylabel("Genus (GTDB)", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.savefig("../output/MAG_top_20_genera.png", dpi=300, bbox_inches="tight")

# Plot the top 20 classes
plt.figure(figsize=(10, 8))
sns.barplot(x=top_20_classes.values, y=top_20_classes.index, color="#96B0CA")

plt.xlabel("Number of MAGs", fontsize=14)
plt.ylabel("Class (based on GTDB)", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.savefig("../output/MAG_top_20_classes.png", dpi=300, bbox_inches="tight")

# Plot the top 20 orders
plt.figure(figsize=(10, 8))
sns.barplot(x=top_20_orders.values, y=top_20_orders.index, color="#96B0CA")

plt.xlabel("Number of MAGs", fontsize=14)
plt.ylabel("Order (based on GTDB)", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12, fontstyle='italic')
plt.bar_label(plt.gca().containers[0], 
                     labels=[f"{x:,}" for x in top_20_orders.values],
                     label_type='edge', 
                     fontsize=10)

plt.savefig("../output/MAG_top_20_orders.png", dpi=300, bbox_inches="tight")
