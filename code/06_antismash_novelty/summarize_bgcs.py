import glob, os, re
import pandas as pd
from Bio import SeqIO

base = "/PATH/TO/virgibacillus_v2"
rows = []
dirs = sorted(d for d in glob.glob(f"{base}/06_antismash/GCA_*") if os.path.isdir(d))
for d in dirs:
    acc = os.path.basename(d)
    for gb in sorted(glob.glob(f"{d}/*.region*.gbk")):
        rec = next(SeqIO.parse(gb, "genbank"))
        for f in rec.features:
            if f.type == "region":
                q = f.qualifiers
                rows.append({
                    "accession": acc, "contig": rec.id, "file": os.path.basename(gb),
                    "region": q.get("region_number", [""])[0],
                    "product": "+".join(q.get("product", [])),
                    "contig_edge": q.get("contig_edge", ["False"])[0] == "True",
                    "length_bp": len(rec.seq)})
b = pd.DataFrame(rows)
b.to_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t", index=False)

qc = pd.read_csv(f"{base}/03_qc/checkm2/quality_report.tsv", sep="\t")[
    ["Name", "Contig_N50", "Total_Contigs", "Genome_Size"]].rename(columns={"Name": "accession"})
g = qc[qc.accession.isin([os.path.basename(d) for d in dirs])].copy()
cnt = b.groupby("accession").agg(n_bgc=("region", "size"), n_edge=("contig_edge", "sum"))
g = g.merge(cnt, on="accession", how="left").fillna({"n_bgc": 0, "n_edge": 0})
g["n_bgc"], g["n_edge"] = g.n_bgc.astype(int), g.n_edge.astype(int)
g.to_csv(f"{base}/06_antismash/bgc_per_genome.tsv", sep="\t", index=False)

print("Genomes analysed:", len(g), "| with >=1 BGC:", int((g.n_bgc > 0).sum()))
print("Total BGC regions:", len(b), "| mean per genome: %.2f | median: %.1f" % (g.n_bgc.mean(), g.n_bgc.median()))
print("BGCs on contig edge: %d (%.1f%%)" % (b.contig_edge.sum(), 100 * b.contig_edge.mean()))
print("\nTop product classes (regions):"); print(b["product"].value_counts().head(15).to_string())
g["N50_bin"] = pd.cut(g.Contig_N50, [0, 20e3, 50e3, 100e3, 500e3, 1e9], labels=["<20k", "20-50k", "50-100k", "100-500k", ">500k"])
x = g.groupby("N50_bin", observed=True).agg(genomes=("accession", "size"), mean_bgc=("n_bgc", "mean"), edge_bgc=("n_edge", "sum"), all_bgc=("n_bgc", "sum"))
x["edge_pct"] = (100 * x.edge_bgc / x.all_bgc).round(1); x["mean_bgc"] = x.mean_bgc.round(2)
print("\nBGC counts and edge fraction by assembly N50:"); print(x.to_string())
print("\nSpearman N50 vs BGC count:", g[["Contig_N50", "n_bgc"]].corr(method="spearman").iloc[0, 1].round(3))
