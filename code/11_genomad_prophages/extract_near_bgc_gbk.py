import re, os
import numpy as np, pandas as pd
from Bio import SeqIO

base = "/PATH/TO/virgibacillus_v2"
out = f"{base}/11_genomad/visualizations_gbk"
os.makedirs(out, exist_ok=True)
norm = lambda s: re.sub(r"\.\d+$", "", str(s))
MARGIN = 10000

d = pd.read_csv(f"{base}/11_genomad/bgc_virus_distance.tsv", sep="\t")
d = d[d["dist"] <= 25000].copy()

v = pd.read_csv(f"{base}/11_genomad/all_viruses.tsv", sep="\t")
v["contig"] = v["seq_name"].str.split("|").str[0].map(norm)
x = v["seq_name"].str.extract(r"\|provirus_(\d+)_(\d+)")
v["pstart"] = pd.to_numeric(x[0]).fillna(0)
v["pend"] = pd.to_numeric(x[1]).fillna(np.inf)

rows = []
for r in d.itertuples():
    g, c = r.genome, r.contig
    gbk = f"{base}/06_antismash/{g}/{g}.gbk"
    rec = next((q for q in SeqIO.parse(gbk, "genbank") if norm(q.id) == c), None)
    if rec is None:
        print("NOT FOUND:", g, c); continue

    s = v[v["contig"] == c].copy()
    s["gap"] = np.maximum(0, np.maximum(s["pstart"] - r.end, r.start - s["pend"]))
    p = s.sort_values("gap").iloc[0]
    ps, pe = int(p["pstart"]), int(min(p["pend"], len(rec)))

    lo = max(0, min(int(r.start), ps) - MARGIN)
    hi = min(len(rec), max(int(r.end), pe) + MARGIN)
    sub = rec[lo:hi]
    sub.annotations["molecule_type"] = "DNA"
    fam = "nofam" if pd.isna(r.fam) else re.sub(r"[^A-Za-z0-9]+", "_", str(r.fam))
    name = f"{fam}__{g}__{c}"
    sub.id = c[:16]; sub.name = c[:16]; sub.description = name
    SeqIO.write(sub, f"{out}/{name}.gbk", "genbank")

    n_cds = sum(f.type == "CDS" for f in sub.features)
    rows.append(dict(file=name, genome=g, contig=c, fam=fam, window=f"{lo}-{hi}",
                     win_len=hi - lo, bgc_in_win=f"{int(r.start)-lo}-{int(r.end)-lo}",
                     provirus_in_win=f"{max(ps,lo)-lo}-{pe-lo}", n_cds=n_cds))
    print(f"{name}: {hi-lo} bp, {n_cds} CDS")

pd.DataFrame(rows).to_csv(f"{out}/windows.tsv", sep="\t", index=False)
