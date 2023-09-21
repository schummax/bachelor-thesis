import pandas as pd
import os

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

path = "/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/snakemake_assembly/binning/metawrap/BIN_REFINEMENT"

total_MAGs = 0

for folder in os.listdir(path):
    if folder.startswith("BGG"):
        for file in os.listdir(os.path.join(path, folder)):
            if file == "metawrap_50_10_bins.stats":
                sample_name = folder.split("-")[0]
                sample_name = sample_name_mapping[sample_name]
                print(f"----------File: {sample_name}----------")
                df = pd.read_csv(os.path.join(path, folder, file), sep="\t")
                rows = df.shape[0]
                print(f"MAGs: {rows}")
                total_MAGs += rows

print(f"Total MAGs: {total_MAGs}")

#############

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

for key, df in dfs.items():
    key = key.replace(".tsv", "")
    # rename with mapping
    key = sample_name_mapping[key]
    print(f"\n----------{key}----------")
    comlete = 0
    incomplete = 0
    for index, row in df.iterrows():
        
        if row['BGC_complete'] == 'Yes':
            comlete += 1
            # total_complete_BGCs += 1
        else:
            incomplete += 1
    print(f"Complete: {comlete}, Incomplete: {incomplete}")