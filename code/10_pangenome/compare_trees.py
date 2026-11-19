import pandas as pd, re, glob
from Bio import Phylo
base = "/PATH/TO/virgibacillus_v2"
res = sorted(glob.glob(f"{base}/10_pangenome/proteins/OrthoFinder/Results_*"))[-1]
reps = pd.read_csv(f"{base}/10_pangenome/b_reps.tsv", sep="\t", dtype=str)
og = [l.strip() for l in open(f"{base}/10_pangenome/outgroups.txt")]
norm = lambda s: re.sub(r"^usr_", "", s)
t_of = Phylo.read(f"{res}/Species_Tree/SpeciesTree_rooted.txt", "newick")
t_gt = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick")
rep_set = set(reps.accession)
print("Tips: OrthoFinder tree %d | GTDB-marker tree %d" % (len(t_of.get_terminals()), len(t_gt.get_terminals())))
print("OrthoFinder tip names (first 3):", [x.name for x in t_of.get_terminals()][:3])

kids = t_of.root.clades
print("\nOrthoFinder root split (tips per side):", [len(k.get_terminals()) for k in kids])
print("Outgroups form one side of the root:", any(set(norm(x.name) for x in k.get_terminals()) == set(og) for k in kids))

def prune_to(t, keep):
    for x in list(t.get_terminals()):
        if norm(x.name) not in keep: t.prune(x)
    return t
t_of = prune_to(t_of, rep_set); t_gt = prune_to(t_gt, rep_set)
names = frozenset(norm(x.name) for x in t_of.get_terminals())
print("Tips shared after pruning to the representatives: OrthoFinder %d | GTDB %d | in common %d" % (
    len(names), len(t_gt.get_terminals()), len(names & frozenset(norm(x.name) for x in t_gt.get_terminals()))))
names = names & frozenset(norm(x.name) for x in t_gt.get_terminals())
ref = sorted(names)[0]
def splits(t):
    out = set()
    for c in t.get_nonterminals():
        s = frozenset(norm(x.name) for x in c.get_terminals()) & names
        if ref in s: s = frozenset(names - s)
        if 1 < len(s) < len(names) - 1: out.add(s)
    return out
A, B = splits(t_of), splits(t_gt)
print("\nSplits: OrthoFinder %d | GTDB-marker %d | shared %d | Robinson-Foulds = %d | normalised RF = %.2f (0 identical, 1 nothing shared)" % (
    len(A), len(B), len(A & B), len(A ^ B), len(A ^ B) / (len(A) + len(B))))
print("\nIs each GTDB species a clade? (OrthoFinder tree / GTDB-marker tree)")
for sp, d in reps.groupby("species"):
    S = frozenset(d.accession) & names
    if len(S) < 2: continue
    S = frozenset(names - S) if ref in S else S
    print("%-30s n=%d  %s / %s" % (sp, len(S), S in A, S in B))
