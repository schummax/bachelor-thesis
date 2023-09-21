import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

# Settings for Seaborn plots
sns.set(font_scale=1.5)  # This multiplies the current font size by 1.5.
sns.set_style("ticks")

# Loading the data
path_to_file = "../../05-results/QUAL_pyDamage_results_BGCcontigs.tsv"
pydamage = pd.read_csv(path_to_file, sep="\t")

# Filtering the data
flt_res = pydamage[
    (pydamage['predicted_accuracy'] >= 0.5) & (pydamage['qvalue'] < 0.05)
][['sample', 'reference', 'nb_reads_aligned'] + [f'CtoT-{i}' for i in range(11)]]

# Melting the data from wide to long format
flt_res = flt_res.melt(id_vars=['sample', 'reference', 'nb_reads_aligned'], 
                       value_vars=[f'CtoT-{i}' for i in range(11)], 
                       var_name='pos', 
                       value_name='freq')

# Extracting the position number and adding 1
flt_res['pos'] = flt_res['pos'].apply(lambda x: int(re.search(r"\d+", x).group()) + 1)

# Filtering rows where nb_reads_aligned is at least 500
flt_res = flt_res[flt_res['nb_reads_aligned'] >= 500]

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
flt_res['sample'] = flt_res['sample'].replace(sample_name_mapping)

# Sorting the data by `sample`
flt_res = flt_res.sort_values(by='sample')

# Calculating the number of unique `reference` values per `sample`
flt_res['n_contigs'] = flt_res.groupby('sample')['reference'].transform('nunique')

# Creating a label for each `sample`
flt_res['label'] = flt_res['sample'] + " (n=" + flt_res['n_contigs'].astype(str) + ")"

# Generating the subplot
g = sns.catplot(
    data=flt_res, 
    x='pos', 
    y='freq', 
    col='label', 
    col_wrap=3, 
    kind='box', 
    height=5, 
    aspect=1.2, 
    dodge=True,
    color='lightgrey',
    showfliers=True, # True = show outliers of the boxplot
)

g.map(sns.stripplot, 'pos', 'freq', jitter=True, color='grey', dodge=True, alpha=0.5)

g.set_axis_labels("position from the 5' end [bp]", "C to T substitution frequency")
g.set_titles("{col_name}")

# Formatting the y-axis to display percentages
for ax in g.axes.flat:
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: '{:.0%}'.format(x)))

plt.tight_layout()
output_path = "../output/pydamage_boxplot.png"
plt.savefig(output_path, format='png', dpi=300, bbox_inches='tight')
plt.show()