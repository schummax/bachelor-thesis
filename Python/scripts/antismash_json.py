import json
import pandas as pd
import os

data_list = []  # A list to hold the data dictionaries

for dir in os.listdir('/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/antismash/'):
    if dir.startswith("BGG_"):
        print(dir)
        tsv_file = pd.read_csv(f'/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/antismash/{dir}/{dir}.tsv', sep='\t')
        with open(f'/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/antismash/{dir}/{dir}.json', 'r') as f:
            print("loading json data")
            json_data = json.load(f)
            for elem in json_data["records"]:
                matching_rows = tsv_file[(tsv_file['Contig_ID'] == elem['id']) & (tsv_file['BGC_complete'] == "Yes")]
                if not matching_rows.empty:
                    cluster_compare_module = elem["modules"].get("antismash.modules.cluster_compare")
                    if cluster_compare_module:
                        res = cluster_compare_module['db_results']['MIBiG']['by_region'].get('1')
                        if res:
                            record_id = elem['id']
                            regions = res['ProtoToRegion_RiQ']['reference_regions'].keys()
                            try:
                                for region in regions:
                                    identity = res['ProtoToRegion_RiQ']['details']['details']['1'][region]['identity']
                                    product = res['ProtoToRegion_RiQ']['reference_regions'][region]['products']
                                    organism = res['ProtoToRegion_RiQ']['reference_regions'][region]['organism']
                                    
                                    # Create a dictionary and add it to the list
                                    data_dict = {
                                        "sample": dir,
                                        "ID": record_id,
                                        "Region": region,
                                        "Identity": identity,
                                        "Product": product,
                                        "Organism": organism
                                    }
                                    data_list.append(data_dict)
                                    print(data_dict)
                            except KeyError as e:
                                print(e)
                                pass

# Convert the list of dictionaries to a DataFrame
df = pd.DataFrame(data_list)

# To save the DataFrame to a CSV file (optional)
df.to_csv('../output/comBGCs_knownclusterblast.tsv', sep="\t", index=False)

# To display the first few rows of the DataFrame (optional)
print(df.head())

# node_count = 0
# for dir in os.listdir('/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/antismash/'):
#     if dir.startswith("BGG_"):
#         print(dir)
#         tsv_file = pd.read_csv(f'/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/antismash/{dir}/{dir}.tsv', sep='\t')
#         with open(f'/mnt/archgen/users/schumacher/bachelorthesis/04-analysis/antismash/{dir}/{dir}.json', 'r') as f:
#             print("loading json data")
#             json_data = json.load(f)
#             for elem in json_data["records"]:
#                 clusterblast_module = elem["modules"].get("antismash.modules.clusterblast")
#                 if clusterblast_module:
#                     knowncluster = clusterblast_module.get("knowncluster")
#                     if knowncluster:
#                         results = knowncluster["results"]
#                         for result in results:
#                             region_number = result["region_number"]
#                             total_hits = result["total_hits"]
#                             if int(total_hits) > 0:
#                                 matching_rows = tsv_file[(tsv_file['Contig_ID'] == elem['id']) & (tsv_file['BGC_complete'] == "Yes")]
#                                 if not matching_rows.empty:
#                                     node_count += 1
#                                     print(f"ID: {elem['id']}")
#                                     print(f"Region {region_number}: Total hits found = {total_hits}")
                                    
#                                     # Create a dictionary and add it to the list
#                                     data_dict = {
#                                         "sample": dir,
#                                         "ID": elem['id'],
#                                         "Region": region_number,
#                                         "Total hits": total_hits,
#                                     }
#                                     data_list.append(data_dict)
#                                     print(data_dict)
                    
# print(f"Total nodes: {node_count}")

# # Convert the list of dictionaries to a DataFrame
# df = pd.DataFrame(data_list)

# # To save the DataFrame to a CSV file (optional)
# df.to_csv('../output/comBGCs_knownclusterblast.tsv', sep="\t", index=False)