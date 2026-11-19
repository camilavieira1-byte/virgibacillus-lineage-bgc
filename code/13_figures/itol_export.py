import pandas as pd, numpy as np, glob, os, re, shutil
from Bio import Phylo
base = "/PATH/TO/virgibacillus_v2"; out = f"{base}/13_figures/itol"; os.makedirs(out, exist_ok=True)
res = sorted(glob.glob(f"{base}/10_pangenome/proteins/OrthoFinder/Results_*"))[-1]
tf = f"{res}/Species_Tree/SpeciesTree_rooted.txt"; shutil.copy(tf, f"{out}/tree_orthofinder.nwk")
t = Phylo.read(tf, "newick"); norm = lambda s: re.sub(r"^usr_", "", s or "")
tips = {norm(x.name): x for x in t.get_terminals()}
G = pd.read_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", dtype=str).fillna("").set_index("accession")
og = {l.strip() for l in open(f"{base}/10_pangenome/outgroups.txt") if l.strip()}
COL = {"Virgibacillus salexigens": "#E8999B", "Virgibacillus salarius": "#F4D6D0", "Virgibacillus pantothenticus": "#F7C99A",
       "Virgibacillus proomii": "#A0467A", "Virgibacillus chiguensis": "#C07BA0", "Virgibacillus dokdonensis": "#EBD0DE", "outgroup": "#D9D9D9"}
def group(a): return "outgroup" if a in og else G.species.get(a, "") or "unclassified"
def label(a):
    s = G.species.get(a, "") or "unclassified"
    return f"{s} | {a}" + (" (outgroup)" if a in og else "")
print("tips in tree:", len(tips), "| outgroups found:", len(og & set(tips)), "| tips without a GTDB species:", sum(1 for a in tips if not G.species.get(a, "") and a not in og))

with open(f"{out}/labels.txt", "w") as f:
    f.write("LABELS\nSEPARATOR TAB\nDATA\n")
    for a, x in tips.items(): f.write(f"{x.name}\t{label(a)}\n")

rows = []
for g, leaves in pd.Series({a: group(a) for a in tips}).groupby(lambda a: group(a)).groups.items() if False else []: pass
grp = {}
for a in tips: grp.setdefault(group(a), []).append(a)
with open(f"{out}/tree_colors.txt", "w") as f:
    f.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
    for g, ac in grp.items():
        col = COL.get(g, "#CCCCCC"); nodes = [tips[a] for a in ac]
        if len(nodes) == 1: f.write(f"{nodes[0].name}\tlabel\t#000000\t{g}\n"); continue
        ca = t.common_ancestor(nodes); inside = {norm(x.name) for x in ca.get_terminals()}
        if inside != set(ac): print("WARNING: %s is not a clade (%d tips in the clade, %d expected)" % (g, len(inside), len(ac)))
        c0, c1 = ca.clades[0], ca.clades[1]
        f.write(f"{c0.get_terminals()[0].name}|{c1.get_terminals()[0].name}\trange\t{col}\t{g}\n")

reg = pd.read_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t"); reg["tok"] = reg["product"].str.split("+")
tk = reg.explode("tok")
classes = ["NRPS", "NRP-metallophore", "RiPP-like", "lanthipeptide-class-ii", "betalactone"]
colors = ["#7A2F5E", "#D9697A", "#F4B183", "#8C6BB1", "#4C78A8"]
pres = {c: set(tk[tk.tok == c].accession) for c in classes}
with open(f"{out}/dataset_bgc_classes.txt", "w") as f:
    f.write("DATASET_BINARY\nSEPARATOR TAB\nDATASET_LABEL\tBGC classes\nCOLOR\t#7A2F5E\n")
    f.write("FIELD_SHAPES\t" + "\t".join(["2"] * len(classes)) + "\nFIELD_LABELS\t" + "\t".join(classes) + "\nFIELD_COLORS\t" + "\t".join(colors) + "\nDATA\n")
    for a, x in tips.items(): f.write(x.name + "\t" + "\t".join("1" if a in pres[c] else "-1" for c in classes) + "\n")
npro = pd.read_csv(f"{base}/11_genomad/mge_per_genome.tsv", sep="\t").set_index("accession").n_prophage
nbgc = reg[reg["product"] != "terpene-precursor"].groupby("accession").size()
for name, ser, col, lab in (("dataset_prophages.txt", npro, "#7A2F5E", "Prophages (geNomad)"), ("dataset_bgc_regions.txt", nbgc, "#D9697A", "BGC regions (no terpene-precursor)")):
    with open(f"{out}/{name}", "w") as f:
        f.write(f"DATASET_SIMPLEBAR\nSEPARATOR TAB\nDATASET_LABEL\t{lab}\nCOLOR\t{col}\nDATA\n")
        for a, x in tips.items(): f.write(f"{x.name}\t{int(ser.get(a, 0))}\n")

low = []
for c in t.get_nonterminals():
    s = c.confidence if c.confidence is not None else c.name
    try: s = float(s)
    except (TypeError, ValueError): continue
    if s < 0.9: low.append((s, len(c.get_terminals())))
print("internal nodes with support < 0.9: %d of %d" % (len(low), len(t.get_nonterminals())), sorted(low)[:8])
print("files written to", out); print(sorted(os.listdir(out)))
