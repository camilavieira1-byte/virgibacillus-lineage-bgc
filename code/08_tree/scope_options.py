import pandas as pd, re
from Bio import Phylo
base = "/PATH/TO/virgibacillus_v2"
s = pd.read_csv(f"{base}/04_gtdbtk_b/gtdb_classification.tsv", sep="\t")
s["tip"] = "usr_" + s.accession
gen = dict(zip(s.tip, s.genus.fillna("NA"))); fam = dict(zip(s.tip, s.family.fillna("NA")))
qc = dict(zip(s.tip, s.pass_qc))
n50 = pd.read_csv(f"{base}/06_antismash/bgc_per_genome.tsv", sep="\t").set_index("accession").Contig_N50.to_dict()

# clone groups (skani, ANI >= 99.9, AF >= 50)
ed = pd.read_csv(f"{base}/05_ani/skani_edges.tsv", sep="\t")
acc = lambda x: re.search(r"(GCA_\d+\.\d)", x).group(1)
ed["a"], ed["b"] = ed.Ref_file.map(acc), ed.Query_file.map(acc)
ed["af"] = ed[["Align_fraction_ref", "Align_fraction_query"]].min(axis=1)
genomes = [acc(l) for l in open(f"{base}/05_ani/genome_list.txt")]
par = {g: g for g in genomes}
def f(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for r in ed[(ed.ANI >= 99.9) & (ed.af >= 50)].itertuples(): par[f(r.a)] = f(r.b)
cg = {g: f(g) for g in genomes}

t = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick")
og = [x.name for x in t.get_terminals() if gen.get(x.name, "").startswith("Paenibacillus")]
t.root_with_outgroup(*[t.find_any(name=n) for n in og])
tips = [x.name for x in t.get_terminals()]

def clade(names):
    ca = t.common_ancestor([t.find_any(name=n) for n in names])
    return [x.name for x in ca.get_terminals()], ca.confidence
def stats(label, names):
    p = [n for n in names if qc.get(n)]
    ids = [n[4:] for n in p]
    print("%-52s all=%3d  passQC=%3d  cloneGroups=%3d  N50>=100kb=%3d" % (
        label, len(names), len(p), len({cg[i] for i in ids if i in cg}),
        sum(1 for i in ids if n50.get(i, 0) >= 100e3)))

vir = [x for x in tips if gen.get(x, "").startswith("Virgibacillus")]
ss = [x for x in tips if gen.get(x) == "Virgibacillus"]
cl, sup = clade(vir)
intr = [x for x in cl if x not in vir]
print("Smallest clade with all Virgibacillus*: %d tips, support %s" % (len(cl), None if sup is None else round(sup, 2)))
print("  of which Virgibacillus*: %d | other genera inside: %d" % (len(cl) - len(intr), len(intr)))
print("  other genera inside:", pd.Series([gen[x] for x in intr]).value_counts().to_dict())
print("  GTDB families among Virgibacillus*:", pd.Series([fam[x] for x in vir]).value_counts().to_dict())
print("  GTDB families among the intruders:", pd.Series([fam[x] for x in intr]).value_counts().to_dict())
print()
stats("A'  clade containing all Virgibacillus*", cl)
stats("B   GTDB Virgibacillus s.s.", ss)
stats("C0  Virgibacillus* only (106, paraphyletic)", vir)
for g in sorted(set(gen[x] for x in vir)):
    stats("    " + g, [x for x in tips if gen.get(x) == g])
