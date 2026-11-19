import re, glob, os
import pandas as pd

base = "/PATH/TO/virgibacillus_v2"
norm = lambda s: re.sub(r"\.\d+$", "", str(s))

# 1. BGC regions: contig + original coordinates + product, from the GBKs
rows = []
for p in glob.glob(f"{base}/07_bigscape/input_fixed/*.gbk"):
    txt = open(p).read()
    head = txt.split("\nFEATURES")[0]
    ver = re.search(r"^VERSION\s+(\S+)", head, re.M).group(1)
    s = re.search(r"Orig\. start\s*::\s*(\d+)", head)
    e = re.search(r"Orig\. end\s*::\s*(\d+)", head)
    prod = re.search(r'/product="([^"]+)"', txt)
    rows.append(dict(file=os.path.basename(p),
                     genome="_".join(os.path.basename(p).split("_")[:2]),
                     contig=norm(ver),
                     start=int(s.group(1)) if s else None,
                     end=int(e.group(1)) if e else None,
                     product=prod.group(1) if prod else "NA"))
bgc = pd.DataFrame(rows)
print("BGC regions:", len(bgc), "| with coordinates:", bgc["start"].notna().sum())

# 2. geNomad tables
def load(path):
    df = pd.read_csv(path, sep="\t")
    df["contig"] = df["seq_name"].str.split("|").str[0].map(norm)
    m = df["seq_name"].str.extract(r"\|provirus_(\d+)_(\d+)")
    df["pstart"] = pd.to_numeric(m[0]); df["pend"] = pd.to_numeric(m[1])
    return df
vir = load(f"{base}/11_genomad/all_viruses.tsv")
pla = load(f"{base}/11_genomad/all_plasmids.tsv")

# 3. Sanity check: do any BGC contigs appear in geNomad at all?
print("BGC contigs found in viruses table :", bgc["contig"].isin(vir["contig"]).sum())
print("BGC contigs found in plasmids table:", bgc["contig"].isin(pla["contig"]).sum())

# 4. Annotate
bgc["on_plasmid"] = bgc["contig"].isin(pla["contig"])

def virus_hit(r):
    sub = vir[vir["contig"] == r["contig"]]
    if sub.empty: return "none"
    if sub["pstart"].isna().any(): return "viral_contig"          # whole contig is viral
    if pd.notna(r["start"]):
        ov = sub[(sub["pstart"] <= r["end"]) & (sub["pend"] >= r["start"])]
        if not ov.empty: return "inside_provirus"
    return "contig_has_provirus_no_overlap"
bgc["virus"] = bgc.apply(virus_hit, axis=1)

bgc["mobility"] = "chromosomal"
bgc.loc[bgc["virus"].isin(["inside_provirus", "viral_contig"]), "mobility"] = "virus"
bgc.loc[bgc["on_plasmid"], "mobility"] = "plasmid"

bgc.to_csv(f"{base}/06_antismash/bgc_mobility_annotated.tsv", sep="\t", index=False)

# 5. Summary
print("\nMobility:\n", bgc["mobility"].value_counts().to_string())
print("\nProvirus relationship:\n", bgc["virus"].value_counts().to_string())
mob = bgc[bgc["mobility"] != "chromosomal"]
if not mob.empty:
    print("\nProduct x mobility:")
    print(pd.crosstab(mob["product"], mob["mobility"]).to_string())
