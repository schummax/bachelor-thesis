import collections
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re

# Settings for Seaborn plots
sns.set(font_scale=1.5)  # This multiplies the current font size by 1.5.
sns.set_style("ticks")

# Loading the pydamage data
path_to_pydamage = "../../05-results/QUAL_pyDamage_results_BGCcontigs.tsv"
pydamage = pd.read_csv(path_to_pydamage, sep="\t")

# Filtering the data
flt_res = pydamage[
    (pydamage['predicted_accuracy'] >= 0.5) & (pydamage['qvalue'] < 0.05)
][['sample', 'reference', 'nb_reads_aligned'] + [f'CtoT-{i}' for i in range(11)]]

# Filtering rows where nb_reads_aligned is at least 500
flt_res = flt_res[flt_res['nb_reads_aligned'] >= 500]

# Loading the MMseqs2 data
path_to_MMseqs = "/mnt/archgen/users/schumacher/bachelorthesis/05-results/ANNO_BGC_contig_taxclassification.tsv"
MMseqs = pd.read_csv(path_to_MMseqs, sep="\t")

# Filtering the data with cotigs that passed pyDamage filtering
filtered_MMseqs = MMseqs[MMseqs['contig'].isin(flt_res['reference'])]

# reset the index
filtered_MMseqs = filtered_MMseqs.reset_index(drop=True)

# Extracting the taxonomic level from the lineage column
filtered_MMseqs["genus"] = filtered_MMseqs["lineage"].apply(lambda x: re.search(r"g_([^;]+)", x).group(1) if re.search(r"g_([^;]+)", x) else np.nan)
filtered_MMseqs["class"] = filtered_MMseqs["lineage"].apply(lambda x: re.search(r"c_([^;]+)", x).group(1) if re.search(r"c_([^;]+)", x) else np.nan)
filtered_MMseqs["order"] = filtered_MMseqs["lineage"].apply(lambda x: re.search(r"o_([^;]+)", x).group(1) if re.search(r"o_([^;]+)", x) else np.nan)

filtered_MMseqs["genus"] = filtered_MMseqs["genus"].replace(np.nan, "unclassified")
filtered_MMseqs["class"] = filtered_MMseqs["class"].replace(np.nan, "unclassified")
filtered_MMseqs["order"] = filtered_MMseqs["order"].replace(np.nan, "unclassified")

# Mapping dictionary
sample_name_mapping = {
    "BGG_BOSW": "0_NA_B",
    "BGG_POSDC": "1_11_P",
    "BGG_KOSDC": "1_16_K",
    "BGG_LOSDA": "2_12_L",
    "BGG_NOSDA": "2_16_N",
    "BGG_ROSDA": "2_19_R",
    "BGG_QOSDA": "3_06_Q",
    "BGG_MOSDA": "3_16_M",
    "BGG_OOSDC": "3_19_O"
}

# Map the old sample names to the new names
filtered_MMseqs['sample'] = filtered_MMseqs['sample'].replace(sample_name_mapping)

# Sorting the data by `sample`
filtered_MMseqs = filtered_MMseqs.sort_values(by='sample')

# print the number of contig per `sample`
n = {}
for sample in filtered_MMseqs["sample"].unique():
    n[sample] = len(filtered_MMseqs[filtered_MMseqs['sample'] == sample])
    print(f"{sample}: {n[sample]}")

# Count the occurrence of each level and select the top 10
top_10_genera = filtered_MMseqs["genus"].value_counts().head(10)
top_10_classes = filtered_MMseqs["class"].value_counts().head(10)
top_10_orders = filtered_MMseqs["order"].value_counts().head(10)

# Plot the top 10 orders
plt.figure(figsize=(10, 8))
sns.barplot(x=top_10_orders.values, y=top_10_orders.index, color="#96B0CA")

plt.xlabel("Number of contigs", fontsize=14)
plt.ylabel("Order", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12, fontstyle='italic')

plt.savefig("../output/MMseqs_top_10_orders.png", dpi=300, bbox_inches="tight")

# list top 8 orders by sample
top_8_orders = filtered_MMseqs.groupby("sample")["order"].value_counts().groupby(level=0).head(8).to_dict()

grouped_orders = collections.defaultdict(dict)
sample_order_count = collections.defaultdict(dict)
for (sample, order_name), count in top_8_orders.items():
    print(f"{sample}: {order_name} ({count})")
    sample_order_count[sample][order_name] = count
# for sample, count in n.items():
#     print(f"{sample}: {count}")
for (sample, order_name), count in top_8_orders.items():
    count = count / n[sample] * 100
    grouped_orders[sample][order_name] = count
    print(f"{sample}: {order_name} ({count:.2f}%)")
    
for sample, orders in grouped_orders.items():
    print(f"{sample}: {orders}")
    
# Setting a constant bar height
bar_height = 0.8

# Plotting the data with the specified modifications
fig, axs = plt.subplots(3, 3, figsize=(22, 15))
fig.subplots_adjust(hspace=0.8, wspace=0.6)

# Determine the maximum number of bars in a subplot to set consistent y-axis limits
max_bars = max([len(orders) for sample, orders in grouped_orders.items()])

for i, (sample, orders) in enumerate(grouped_orders.items()):
    if i < 9:  # We have only 9 subplots
        ax = axs.flat[i]
        order_names = list(orders.keys())
        counts = list(orders.values())
        n_contigs = n[sample]
        absolute_counts = [sample_order_count[sample][order_name] for order_name in order_names]
        sns.barplot(x=counts, y=order_names, color="#96B0CA", ax=ax, orient="h", dodge=False, ci=None)
        ax.set_title(f"{sample} (n={n_contigs})", fontsize=14, pad=20)
        ax.set_xlabel("Rel. abundance", fontsize=14)
        ax.set_ylabel("Order", fontsize=14)
        ax.tick_params(axis="x", labelsize=12, )
        ax.tick_params(axis="y", labelsize=12, width=bar_height)
        for label in ax.get_yticklabels():
            if label != "unclassified":
                label.set_fontstyle('italic')
        
        # Adjusted bar labels to include absolute counts along with percentages
        ax.bar_label(ax.containers[0], 
                     labels=[f"{count:.2f}% ({abs_count})" for abs_count, count in zip(absolute_counts, counts)], 
                     label_type='edge', 
                     fontsize=10)
        ax.set_ylim(-1, max_bars)
        ax.set_xlim(0, 100)

plt.savefig("../output/MMseqs_top_8_orders_by_sample.png", dpi=300, bbox_inches="tight")