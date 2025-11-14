import pandas as pd
import networkx as nx

call_samples = pd.read_csv("plink.imiss", sep="\s+")
related_sample = pd.read_csv("pihat_min0.2.genome", sep="\s+")

# Prepare sample groups and F_miss dict for samples
G = nx.from_pandas_edgelist(related_sample, 'IID1', 'IID2')
groups = list(nx.connected_components(G))
F_dict = dict(call_samples[["IID", "F_MISS"]].values)

# Make sure sample contains in miss file
filter_groups = []
for group in groups:
    correct_group = []
    for sample in group:
        if F_dict.get(sample):
            correct_group.append(sample)
    filter_groups.append(correct_group)

# Find samples should be drop (max F_MISS in group)
drop_samples_list = []
for group in filter_groups:
    scores = list(map(lambda sample: F_dict[sample], group))
    remain_index = scores.index(min(scores))
    drop_samples = group[:remain_index] + group[remain_index+1:]
    drop_samples_list.append(drop_samples)

# Summarize sample groups and samples should be drop
group_summary = pd.DataFrame(data={"groups": filter_groups, "drop_samples": drop_samples_list})
group_summary = group_summary[group_summary["groups"].map(lambda x: len(x) != 1)]

## Create list of dropped samples
drop_samples = group_summary.explode(column='drop_samples')\
                            .merge(call_samples, how="left", left_on="drop_samples", right_on="IID")\
                            .loc[:,["FID", "drop_samples"]]

drop_samples.to_csv("0.2_low_call_rate_pihat.txt", header=None, index=False, sep="\t")
