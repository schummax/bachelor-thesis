import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
dnamollen = pd.read_csv("/mnt/archgen/users/schumacher/bachelorthesis/05-results/QUAL_DNA_molecule_length_dist.tsv", sep="\t")

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

# Update the sample names in the DataFrame
dnamollen['sample'] = dnamollen['sample'].map(sample_name_mapping)

# Calculate fraction
dnamollen['total'] = dnamollen.groupby('sample')['count'].transform('sum')
dnamollen['frac'] = (dnamollen['count'] / dnamollen['total']) * 100

# Specify the order of samples for the legend
sample_order = ['0_NA_B', '1_11_P', '1_16_K', '2_12_L', '2_16_N', '2_19_R', '3_06_Q', '3_16_M', '3_19_O']

sns.set(font_scale=1.1)

# Plot
fig, ax = plt.subplots(figsize=(10, 6))  # Get the Axes object
sns.lineplot(data=dnamollen[dnamollen['length'] <= 290], x='length', y='frac', hue='sample', hue_order=sample_order, errorbar=None, ax=ax)
plt.xlabel("Length [bp]")
plt.ylabel("Fraction of DNA molecules [%]")
leg = plt.legend(title='Sample', loc='upper left')  # Get the legend object
leg.get_frame().set_facecolor('white')
plt.xlim(30, 300)
plt.tight_layout()

ax.set_facecolor('white')

for spine in ax.spines.values():
    spine.set_edgecolor('black')

plt.savefig('length_dist.jpg', dpi=600) 

dnamollen_mean = {}
sample_count = {}
sample_length = {}
# iterate over the rows of the DataFrame and multiply the count by the length
for sample in sample_order:
    print(sample)
    dnamollen_mean[sample] = 0
    sample_count[sample] = 0
    sample_length[sample] = 0
    for index, row in dnamollen.iterrows():
        if row['sample'] == sample and row['length'] != 512: # exclude the 512 bp molecules as they are not real (possible artefact)
            sample_count[sample] += row['count']
            sample_length[sample] += row['count'] * row['length']
        else:
            continue
    #calculate the mean
    dnamollen_mean[sample] = sample_length[sample] / sample_count[sample]

    print(dnamollen_mean[sample])
