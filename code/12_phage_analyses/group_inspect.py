import pandas as pd, os, glob, csv, re, collections
base = "/PATH/TO/virgibacillus_v2"; ph = f"{base}/12_phages"
P = pd.read_csv(f"{ph}/prophages.tsv", sep="\t")
C = pd.read_csv(f"{ph}/clusters_phabox_aai.tsv", sep="\t", header=None, names=["pid", "c"])
G = pd.read_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", dtype=str).fillna("").set_index("accession")
reps = set(pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str).accession)
struct = re.compile(r"capsid|terminase|portal|tail|baseplate|head|tape measure", re.I)
for gname in os.environ.get("GROUPS", "group_1 group_3").split():
    d = P.merge(C, on="pid"); d = d[(d.c == gname) & d.accession.isin(reps)].copy()
    cnt, tot, st = collections.Counter(), collections.Counter(), collections.Counter()
    for acc, sub in d.groupby("accession"):
        f = glob.glob(f"{base}/11_genomad/usr_{acc}/usr_{acc}_find_proviruses/usr_{acc}_provirus_genes.tsv")[0]; want = set(sub.pid)
        for row in csv.reader(open(f), delimiter="\t"):
            if row[0] == "gene": continue
            pid = acc + "__" + row[0].rsplit("_", 1)[0].replace("|", "_")
            if pid not in want: continue
            tot[pid] += 1; desc = row[-1]
            if desc not in ("NA", ""):
                cnt[desc] += 1
                if struct.search(desc): st[pid] += 1
    d["species"] = d.accession.map(lambda a: G.species.get(a, "")); d["genes"] = d.pid.map(tot); d["structural_genes"] = d.pid.map(st).fillna(0).astype(int)
    pd.set_option("display.width", 200)
    print("\n================ %s: %d prophages in representatives" % (gname, len(d)))
    print(d[["accession", "species", "length", "genes", "structural_genes", "n_hallmarks", "in_seq_edge"]].sort_values("length").to_string(index=False))
    print("most common annotations:"); [print("  %3d  %s" % (n, s[:90])) for s, n in cnt.most_common(14)]
