import os
import pandas as pd

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

def get_subdirs(directory):
    """
    Get a list of all sub-directories within a given directory.
    """
    dirs = [os.path.join(directory, d) for d in os.listdir(directory) if os.path.isdir(os.path.join(directory, d))]
    print(f"Found {len(dirs)} sub-directories.")
    
    dirs = [dir for dir in dirs if os.path.basename(dir).startswith('BGG_')]
    print(f"Found {len(dirs)} sub-directories starting with 'BGG_'.")
    
    return dirs

def get_tsv_file(directory):
    """
    Get a list of all TSV files within a given directory.
    """
    files = os.listdir(directory)
    for file in files:
        if file.startswith('BGG_') and file.endswith('.tsv'):
            print(f"Found TSV file: {file}")
            return os.path.join(directory, file)


parent_dir = "/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/antismash"

dfs = {}  # Dictionary to store file path as key and dataframe as value

for sub_directory in get_subdirs(parent_dir):
    tsv_file = get_tsv_file(sub_directory)
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        dfs[os.path.basename(tsv_file)] = df
    except Exception as e:
        print(f"Error reading {tsv_file}. Error: {e}")

table = pd.DataFrame(columns=["Sample", "Contig_ID", "Product_class", "Taxonomic_rank", "Taxonomic_classification", "knownclusterblast"])
taxonomix_results = pd.read_csv("/mnt/archgen/users/schumacher/bachelorthesis/05-results/ANNO_BGC_contig_taxclassification.tsv", sep="\t")

# for key, df in dfs.items():
#     key = key.replace(".tsv", "")
#     key = sample_name_mapping[key]
#     print(f"\n----------{key}----------")
#     for index, row in df.iterrows():
#         if row["BGC_complete"] == "Yes":
#             table["Sample"] = key
#             print(key)
#             table["Contig_ID"] = row["Contig_ID"]
#             print(row["Contig_ID"])
#             table["Product_class"] = row["Product_class"]
#             print(row["Product_class"])
#             table["Taxonomic_classification"] = taxonomix_results[taxonomix_results["contig"] == row["Contig_ID"]]["lineage"]
#             print(taxonomix_results[taxonomix_results["contig"] == row["Contig_ID"]]["lineage"])

# table.to_csv("/mnt/archgen/users/schumacher/bachelorthesis/Python/output/comp_BGC_taxo.tsv", sep="\t", index=False)

# Assuming table and taxonomix_results are already initialized

knownclusterblast_results = pd.read_csv("/mnt/archgen/users/schumacher/bachelorthesis/Python/output/comBGCs_knownclusterblast.tsv", sep="\t")

def knownclusterblast(id):
    if id in knownclusterblast_results["ID"].values:
        print(f"ID {id} FOUND ")
        # sum all rows with the same ID
        count = 0
        for index, row in knownclusterblast_results.iterrows():
            if row["ID"] == id:
                count += row["Total hits"]
        return count
    else:
        print(f"ID {id} not found in knownclusterblast_results")
        return 0

rows_list = []  # Use a list to collect rows and append them to the dataframe at once

for key, df in dfs.items():
    key = key.replace(".tsv", "")
    key = sample_name_mapping[key]
    
    for _, row in df.iterrows():
        if row["BGC_complete"] == "Yes":
            lineage_series = taxonomix_results[taxonomix_results["contig"] == row["Contig_ID"]]["lineage"]
            rank = taxonomix_results[taxonomix_results["contig"] == row["Contig_ID"]]["NCBIrank"]
            
            if not lineage_series.empty:
                lineage = lineage_series.iloc[0]
            else:
                lineage = None  # Or some default value
                
            if not rank.empty:
                rank = rank.iloc[0]
            else:
                rank = None
                
            # Create a new row as a dictionary and append to rows_list
            rows_list.append({
                "Sample": key,
                "Contig_ID": row["Contig_ID"],
                "Product_class": row["Product_class"],
                "Taxonomic_rank": rank,
                "Taxonomic_classification": lineage,
                "knownclusterblast": knownclusterblast(row["Contig_ID"])
            })

# Convert rows_list to a dataframe and append it to table
new_df = pd.DataFrame(rows_list)
table = table._append(new_df, ignore_index=True)

# Save to CSV
table.to_csv("/mnt/archgen/users/schumacher/bachelorthesis/Python/output/comp_BGC_taxo.tsv", sep="\t", index=False)


print(new_df["Sample"].value_counts())
print(new_df["Product_class"].value_counts())
print(new_df["knownclusterblast"].value_counts())

for sample in new_df["Sample"].unique():
    print(f"\n{sample}")
    print(new_df[new_df["Sample"] == sample]["Product_class"].value_counts())

top_knowclusterblast = new_df.sort_values(by="knownclusterblast", ascending=False).head(10)
print(top_knowclusterblast)


