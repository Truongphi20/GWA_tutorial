import pandas as pd

okgp_file = "1kG_MDS7.bim"
hapmap_file = "HapMap-adj.bim"

col_names = ["chr", "rsID", "GD", "pos", "ref", "alt"]

## Create common table
bim1_df = pd.read_csv(okgp_file, sep="\s+", names=col_names)
bim2_df = pd.read_csv(hapmap_file, sep="\s+", names=col_names)
merge_table = bim1_df.merge(bim2_df, how="inner", on="rsID", suffixes=["_okgp", "_hapmap"])

merge_table["okgp_concat"] = (merge_table[f"ref_okgp"] + merge_table[f"alt_okgp"]).map(lambda x: "".join(sorted(x)))
merge_table["hapmap_concat"] = (merge_table[f"ref_hapmap"] + merge_table[f"alt_hapmap"]).map(lambda x: "".join(sorted(x)))

## Check strand flipping
complement_dict = {
    "AC": "GT",
    "AG": "CT",
    "AT": "AT",
    "CG": "CG",
    "CT": "AG",
    "GT": "AC"
}

def check_flipping(snp1, snp2):
    if snp1 == snp2:
        return False
    
    if complement_dict[snp1] == snp2:
        return True
    
    return False

merge_table["check_flipping"] = merge_table[["okgp_concat","hapmap_concat"]]\
                                        .apply(lambda x: check_flipping(x["okgp_concat"], x["hapmap_concat"]), axis="columns")

flipping_table = merge_table[merge_table["check_flipping"]]
flipping_table[["rsID"]].to_csv("flipping_snps.txt", index=None, header=None)

## Get a list of dropped SNPs
drop_table = merge_table[~merge_table["check_flipping"]]
drop_table = merge_table[~merge_table["check_flipping"]].copy()
drop_table["check_dropping"] = drop_table["okgp_concat"] != drop_table["hapmap_concat"]
drop_table = drop_table[drop_table["check_dropping"]]
drop_table["rsID"].to_csv("dropping_snps.txt", index=None, header=None)