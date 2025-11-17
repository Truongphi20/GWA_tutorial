import pandas as pd 

panel_file = "20100804.ALL.panel"
mds_file = "MDS_merge2.mds"
hapmap_file = "HapMap_MDS.fam"

## Read data
race_dict = {
    "JPT": "ASN",
    "ASW": "AFR",
    "CEU": "EUR",
    "CHB": "ASN",
    "CHD": "ASN",
    "YRI": "AFR",
    "LWK": "AFR",
    "TSI": "EUR",
    "MXL": "AMR",
    "GBR": "EUR",
    "FIN": "EUR",
    "CHS": "ASN",
    "PUR": "AMR"
}

okgp_collect = {"IID": [], "race": []}
with open(panel_file, "r") as f:
    for line in f:
        okgp_collect["IID"].append(line.split("\t")[0])
        okgp_collect["race"].append(race_dict[line.split("\t")[1]])
okgp_data = pd.DataFrame(okgp_collect)

mds_data = pd.read_csv(mds_file, sep='\s+')
hapmap_data = pd.read_csv(hapmap_file, sep='\s+', header=None).loc[:,[0,1]].rename(columns={0:"FID", 1:"IID"})

## Get FID of smaples in the panel from mds file 
fill_fid_data = okgp_data.merge(mds_data[["IID", "FID"]], how="left", on="IID")

## Concatenate with OWN data
hapmap_data["race"] = "OWN"
concat_data = pd.concat([fill_fid_data, hapmap_data], axis="index").loc[:,["FID", "IID", "race"]]

## Exporting
concat_data.to_csv("racefile.txt", sep=" ", index=False)