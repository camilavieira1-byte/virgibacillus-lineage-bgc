import pandas as pd, re
base = "/PATH/TO/virgibacillus_v2"; ph = f"{base}/12_phages"
P = pd.read_csv(f"{ph}/prophages.tsv", sep="\t")
C = pd.read_csv(f"{ph}/clusters_phabox_ani.tsv", sep="\t", header=None, names=["pid", "c"])
G = pd.read_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", dtype=str).fillna("")[["accession", "genus", "species", "Contig_N50"]]
d = P.merge(C, on="pid").merge(G, on="accession", how="left")
sp = d[d.species != ""].groupby("c").species.nunique()
x = d[d.c.isin(sp[sp > 1].index)]
pd.set_option("display.width", 220)
print(x[["c", "accession", "species", "length", "in_seq_edge", "virus_score", "n_hallmarks", "Contig_N50"]].to_string(index=False))
H = pd.read_csv(f"{base}/05_ani/skani_edges.tsv", sep="\t")
acc = lambda s: re.search(r"(GCA_\d+\.\d)", s).group(1)
H["a"], H["b"] = H.Ref_file.map(acc), H.Query_file.map(acc)
s = set(x.accession)
print("\nhost ANI between these genomes:")
print(H[H.a.isin(s) & H.b.isin(s)][["a", "b", "ANI"]].to_string(index=False))
