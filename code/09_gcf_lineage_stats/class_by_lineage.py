import pandas as pd, numpy as np
base = "/PATH/TO/virgibacillus_v2"
cg = pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str)
cg = cg[cg.in_A == "True"].set_index("accession")
reg = pd.read_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t")
reg = reg[reg.accession.isin(cg.index)]
tok = reg.assign(tok=reg["product"].str.split("+")).explode("tok")
pres = pd.crosstab(tok.accession, tok.tok).clip(upper=1).reindex(cg.index).fillna(0)
cg["lineage"] = np.where(cg.in_B == "True", "B: " + cg.species.str.replace("Virgibacillus ", ""), cg.genus)
for k in ["NRPS", "RiPP-like", "NRP-metallophore", "lanthipeptide-class-ii"]:
    if k not in pres.columns: continue
    s = pres[k].groupby(cg.lineage).agg(carriers="sum", genomes="size")
    s = s[s.carriers > 0].sort_values("carriers", ascending=False)
    print("\n%s: carried in %d lineages (of %d); in B species: %d of 6" % (k, len(s), cg.lineage.nunique(), int(s.index.str.startswith("B:").sum())))
    print(s.astype(int).to_string())
