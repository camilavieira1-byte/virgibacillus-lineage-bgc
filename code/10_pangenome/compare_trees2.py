import pandas as pd, re, glob
from Bio import Phylo
base = "/PATH/TO/virgibacillus_v2"
res = sorted(glob.glob(f"{base}/10_pangenome/proteins/OrthoFinder/Results_*"))[-1]
reps = pd.read_csv(f"{base}/10_pangenome/b_reps.tsv", sep="\t", dtype=str); keep = set(reps.accession)
norm = lambda s: re.sub(r"^usr_", "", s or "")
def load(p):
    t = Phylo.read(p, "newick")
    for x in list(t.get_terminals()):
        if norm(x.name) not in keep: t.prune(x)
    return t
tf = load(f"{res}/Species_Tree/SpeciesTree_rooted.txt"); tg = load(f"{base}/08_tree/virgi_iqtree.treefile")
names = frozenset(norm(x.name) for x in tf.get_terminals()) & frozenset(norm(x.name) for x in tg.get_terminals())

def splits(t, names):
    ref = sorted(names)[0]; d = {}
    for c in t.get_nonterminals():
        if c is t.root: continue
        s = frozenset(norm(x.name) for x in c.get_terminals()) & names
        if ref in s: s = frozenset(names - s)
        if 1 < len(s) < len(names) - 1:
            lab = c.name if c.name else (str(c.confidence) if c.confidence is not None else "")
            nums = [float(v) for v in re.findall(r"\d+(?:\.\d+)?", lab)]
            d[s] = nums[-1] if nums else None
    return d
SF, SG = splits(tf, names), splits(tg, names)
print("Splits: OrthoFinder %d | GTDB-marker %d | shared %d" % (len(SF), len(SG), len(set(SF) & set(SG))))
print("Support labels, OrthoFinder tree (first 8):", list(SF.values())[:8])
print("Support labels, GTDB-marker tree (first 8, last number = UFBoot):", list(SG.values())[:8])
for thr in (70, 90, 95):
    s = [k for k, v in SG.items() if v is not None and v >= thr]
    print("GTDB-marker splits with UFBoot >= %d: %d | present in the OrthoFinder tree: %d" % (thr, len(s), sum(k in SF for k in s)))
vals = [v for v in SF.values() if v is not None]
if vals:
    sc = 100 if max(vals) <= 1.0 else 1
    for thr in (70, 90, 95):
        s = [k for k, v in SF.items() if v is not None and v * sc >= thr]
        print("OrthoFinder splits with support >= %d (%s scale): %d | present in the GTDB-marker tree: %d" % (thr, "0-1 x100" if sc == 100 else "0-100", len(s), sum(k in SG for k in s)))

sp = reps.groupby("species").accession.first(); sp_names = frozenset(sp.values); lab = {a: s.replace("Virgibacillus ", "V.") for s, a in sp.items()}
print("\nRelationships among the 6 species (one tip per species; species are clades in both trees):")
for nm, t in (("OrthoFinder", tf), ("GTDB-marker", tg)):
    sp_splits = splits(t, sp_names)
    print(" ", nm, ":", sorted(["+".join(sorted(lab[a] for a in s)) for s in sp_splits]))
