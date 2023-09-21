import pandas as pd
import os


def extract_statistics_from_data(data, sample_name_mapping):
    """
    Extracts statistics from the provided data frame.
    """
    # Replace sample names
    data["sample"] = data["sample"].replace(sample_name_mapping)

    # Get sample name
    sample = data["sample"].iloc[0]
    t_sample = sample.replace("_", r"\_")  # replace _ with \_

    # Count pass.MIMAG_high
    MIMAG_high = data[data["pass.MIMAG_high"] == True].shape[0]
    t_MIMAG_high = r"\numprint{" + str(MIMAG_high) + "}"

    # Count pass.MIMAG_medium
    MIMAG_medium = data[data["pass.MIMAG_medium"] == True].shape[0]
    t_MIMAG_medium = r"\numprint{" + str(MIMAG_medium) + "}"

    # Count pass.GUNC
    GUNC = data[data["pass.GUNC"] == True].shape[0]
    t_GUNC = r"\numprint{" + str(GUNC) + "}"

    # Max checkM.completeness
    max_completeness = data["checkM.completeness"].max()
    t_max_completeness = r"\numprint{" + str(max_completeness) + "}"

    # Mean checkM.completeness
    mean_completeness = data["checkM.completeness"].mean().round(2)
    t_mean_completeness = r"\numprint{" + str(mean_completeness) + "}"

    # Format line for LaTeX
    line = (
        r" & ".join(
            [
                t_sample,
                t_MIMAG_high,
                t_MIMAG_medium,
                t_GUNC,
                t_max_completeness,
                t_mean_completeness,
            ]
        )
        + r" \\"
    )

    return line


# Define the sample_name_mapping
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

# Define the LaTeX table formatting
table_top = r"""\begin{table}[H]
\centering
\caption{MAG statistics for all samples.}
\begin{tabular}{ccccccc}
\hline
\textbf{Sample} & \textbf{pass.MIMAG_high} & \textbf{pass.MIMAG_medium} & \textbf{pass.GUNC} & \textbf{Max completeness} & \textbf{Mean completeness} \\
\hline
"""

table_bottom = r"""
\hline
\label{tab:MAG_statistics}
\end{tabular}
\end{table}
"""

path = "/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/snakemake_assembly/stats/bin_quality"

lines = []

for file in os.listdir(path):
    if file.endswith(".tsv"):
        print(f"----------File: {file}----------")
        data = pd.read_csv(os.path.join(path, file), sep="\t", header=0)
        lines.append(extract_statistics_from_data(data, sample_name_mapping))

# Sort lines based on sample_name_mapping order
sorted_lines = sorted(
    lines,
    key=lambda x: list(sample_name_mapping.values()).index(
        x.split(" & ")[0].replace(r"\_", "_")
    ),
)

# Generate the LaTeX table
table_content = "\n".join(sorted_lines)
latex_table = table_top + table_content + table_bottom
print(latex_table)

with open("../output/latex_bin_quality.txt", "w") as file:
    file.write(latex_table)
