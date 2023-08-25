import glob
import matplotlib.pyplot as plt
import pandas as pd

# Define the folder path
folder_path = "/mnt/archgen/users/schumacher/bachelorthesis/tmp/mmseqs_pydamage"

# Mapping dictionary
sample_name_mapping = {
    "BGG_BOSW": "0_NA_B",
    "BGG_POSDC": "1_11_P",
    "BGG_KOSDC": "1_16_K",
    "BGG_LOSDA": "2_12_L",
    "BGG_ROSDA": "2_19_R",
    "BGG_QOSDA": "3_06_Q",
    "BGG_MOSDA": "3_16_M",
}

plt.figure(figsize=(10, 6))

# Iterate through the sample_name_mapping dictionary to maintain order
for sample_name_key, sample_label in sample_name_mapping.items():
    file_path = f"{folder_path}/{sample_name_key}.pydamage.tsv"

    print(f"Processing {file_path}")

    # Read the data
    data = pd.read_csv(file_path, sep=",")

    # Extract the C-to-T transition columns
    ctot_columns = [f"CtoT-{i}" for i in range(20)]

    # Rename the columns to remove the "CtoT-" prefix
    renamed_columns = {col: col.split("-")[-1] for col in ctot_columns}
    data.rename(columns=renamed_columns, inplace=True)

    # Extract the renamed C-to-T transitions
    ctot_data = data[renamed_columns.values()].mean()

    # Plot the C-to-T transitions
    plt.plot(ctot_data, linestyle="-", label=sample_label)

plt.xlabel("Base position from 5'")
plt.ylabel("Substitution frequency")
plt.legend()
plt.xlim(0, 19)
plt.tight_layout()
plt.grid(True)
plt.savefig("../output/pydamage_plot.png", dpi=300)
