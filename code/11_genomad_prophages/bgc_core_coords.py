import re, glob, os
import pandas as pd

base = "/PATH/TO/virgibacillus_v2"
rows = []
for p in glob.glob(f"{base}/07_bigscape/input_fixed/*.gbk"):
    txt = open(p).read()
    head = txt.split("\nFEATURES")[0]
    s = re.search(r"Orig\. start\s*::\s*(\d+)", head)
    if not s: continue
    off = int(s.group(1))
    lo, hi, n = [], [], 0
    for m in re.finditer(r"^\s{5}proto_core\s+(\S+)", txt, re.M):
        nums = [int(x) for x in re.findall(r"\d+", m.group(1))]
        if nums:
            lo.append(min(nums)); hi.append(max(nums)); n += 1
    if n:
        rows.append(dict(stem=os.path.basename(p)[:-4],
                         core_start=off + min(lo) - 1, core_end=off + max(hi), n_cores=n))
c = pd.DataFrame(rows)
c.to_csv(f"{base}/06_antismash/bgc_core_coords.tsv", sep="\t", index=False)
print("regions with proto_core:", len(c), "| with >1 core (hybrids):", (c.n_cores > 1).sum())

a = pd.read_csv(f"{base}/06_antismash/bgc_mobility_prevalence.tsv", sep="\t")
a["stem"] = a["file"].str.replace(r"\.gbk$", "", regex=True)
a = a.merge(c, on="stem", how="left")
has = a["core_start"].notna()
print("BGCs with core coords:", has.sum(), "of", len(a))
a["region_start"], a["region_end"] = a["start"], a["end"]
a.loc[has, "start"] = a.loc[has, "core_start"]
a.loc[has, "end"] = a.loc[has, "core_end"]
print("median length, region vs core (bp):",
      int((a["region_end"] - a["region_start"]).median()),
      int((a["end"] - a["start"]).median()))
a.to_csv(f"{base}/06_antismash/bgc_mobility_prevalence_core.tsv", sep="\t", index=False)
