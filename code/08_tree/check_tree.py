import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo
base = "/PATH/TO/virgibacillus_v2"
s = pd.read_csv(f"{base}/04_gtdbtk_b/gtdb_classification.tsv", sep="\t")
s["tip"] = "usr_" + s.accession
gen = dict(zip(s.tip, s.genus.fillna("NA")))
sp = dict(zip(s.tip, s.species.fillna("")))

t = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick")
tips = [x.name for x in t.get_terminals()]
print("Tips in tree:", len(tips), "| genomes in GTDB table:", len(s))
miss = sorted(set(s.tip) - set(tips))
print("In GTDB table but NOT in tree:", len(miss), miss[:10])

# root on the most distant genera available
for og in ("Paenibacillus_AQ", "Metabacillus", "Terribacillus", "Sediminibacillus", "Oceanobacillus"):
    names = [x for x in tips if gen.get(x, "").startswith(og)]
    if names:
        t.root_with_outgroup(*[t.find_any(name=n) for n in names])
        print("Rooted with:", og, "(", len(names), "genomes )"); break
else:
    print("No outgroup genus found - tree left as is")

def mono(names):
    names = [n for n in names if n in tips]
    if len(names) < 2: return ("single", None, [])
    ca = t.common_ancestor([t.find_any(name=n) for n in names])
    inside = [x.name for x in ca.get_terminals()]
    intruders = [x for x in inside if x not in names]
    sup = ca.confidence
    return ("MONOPHYLETIC" if not intruders else "NOT monophyletic", sup, intruders)

groups = {
 "Virgibacillus sensu lato (all Virgibacillus*)": [x for x in tips if gen.get(x, "").startswith("Virgibacillus")],
 "Virgibacillus s.s. (GTDB)": [x for x in tips if gen.get(x) == "Virgibacillus"],
 "Lentibacillus* (all)": [x for x in tips if gen.get(x, "").startswith("Lentibacillus")],
}
for gname in sorted(set(v for v in gen.values() if v.startswith("Virgibacillus_") or v.startswith("Lentibacillus_"))):
    groups[gname] = [x for x in tips if gen.get(x) == gname]
print("\n%-48s %5s  %-18s %-8s %s" % ("group", "n", "result", "support", "intruders"))
rows = []
for k, v in groups.items():
    r, sup, intr = mono(v)
    print("%-48s %5d  %-18s %-8s %s" % (k, len(v), r, "" if sup is None else round(sup, 2),
          ", ".join(f"{gen.get(i,'?')}:{i[-13:]}" for i in intr[:6]) + (" ..." if len(intr) > 6 else "")))
    rows.append({"group": k, "n": len(v), "result": r, "support": sup, "n_intruders": len(intr), "intruders": ";".join(intr)})
pd.DataFrame(rows).to_csv(f"{base}/08_tree/monophyly_check.tsv", sep="\t", index=False)

# figure
cmap = plt.get_cmap("tab20")
ug = sorted(set(gen.get(x, "NA") for x in tips))
col = {g: matplotlib.colors.to_hex(cmap(i % 20)) for i, g in enumerate(ug)}
fig = plt.figure(figsize=(14, 32)); ax = fig.add_subplot(1, 1, 1)
Phylo.draw(t, axes=ax, do_show=False, show_confidence=False,
           label_func=lambda c: f"{gen.get(c.name,'')} | {c.name[4:]}" if c.is_terminal() else "",
           label_colors=lambda n: col.get(gen.get(n, "NA"), "black"))
ax.set_title("IQ-TREE of 120 bacterial marker genes (GTDB-Tk), tips colored by GTDB genus")
plt.tight_layout(); plt.savefig(f"{base}/08_tree/tree_check.png", dpi=150)
print("\nWrote tree_check.png and monophyly_check.tsv")
