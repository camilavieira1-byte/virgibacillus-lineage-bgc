import pandas as pd, re, os
from Bio import Phylo
base = "/PATH/TO/virgibacillus_v2"; out = f"{base}/13_figures/itol_fig1"; os.makedirs(out, exist_ok=True)
t = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick")
gt = pd.read_csv(f"{base}/04_gtdbtk_b/gtdb_classification.tsv", sep="\t").drop_duplicates("user_genome")
gen = dict(zip(gt.user_genome, gt.genus.fillna(""))); sp = dict(zip(gt.user_genome, gt.species.fillna("")))
tips = {x.name: x for x in t.get_terminals()}
og = [tips[n] for n in gt.user_genome if gen.get(n, "").startswith("Paenibacillus") and n in tips]
t.root_with_outgroup(*og); tips = {x.name: x for x in t.get_terminals()}
def supp(c):
    lab = c.name or (str(c.confidence) if c.confidence is not None else ""); nums = re.findall(r"\d+(?:\.\d+)?", lab)
    return (float(nums[0]), float(nums[-1])) if len(nums) >= 2 else (None, float(nums[0]) if nums else None)
vir = [n for n in tips if gen.get(n, "").startswith("Virgibacillus")]
sets = {"GTDB Virgibacillus": [n for n in tips if gen.get(n) == "Virgibacillus"], "all Virgibacillus-genus genomes": vir}
for g in sorted({gen[n] for n in tips if gen.get(n, "").startswith(("Virgibacillus_", "Lentibacillus"))}): sets[g] = [n for n in tips if gen.get(n) == g]
print("%-34s %5s %9s %-6s %s" % ("group", "tips", "node tips", "mono", "SH-aLRT / UFBoot"))
for k, names in sets.items():
    if len(names) < 2: continue
    c = t.common_ancestor([tips[n] for n in names]); n_in = len(c.get_terminals()); a, b = supp(c)
    print("%-34s %5d %9d %-6s %s / %s" % (k, len(names), n_in, n_in == len(names), a, b))
for c in t.get_nonterminals(): a, b = supp(c); c.confidence = b; c.name = None
Phylo.write(t, f"{out}/tree_marker138.nwk", "newick")
vs = sorted({g for g in gen.values() if g.startswith("Virgibacillus_")}); ln = sorted({g for g in gen.values() if g.startswith("Lentibacillus")})
pv = ["#C07BA0", "#E8999B", "#F4B183", "#F7C99A", "#8C6BB1", "#A0467A", "#EBD0DE", "#D9697A", "#B7A4C9", "#F4D6D0", "#9E6F8E"]; pl = ["#4C78A8", "#8FA1B3", "#3B4252"]
col = {"Virgibacillus": "#7A2F5E"}; col.update({g: pv[i % len(pv)] for i, g in enumerate(vs)}); col.update({g: pl[i % len(pl)] for i, g in enumerate(ln)})
shown = [g for g in col if any(gen.get(n) == g for n in tips)]
with open(f"{out}/dataset_genus.txt", "w") as f:
    f.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tGTDB genus\nCOLOR\t#7A2F5E\n")
    f.write("LEGEND_TITLE\tGTDB r232 genus\nLEGEND_SCALE\t1\nLEGEND_POSITION_X\t60\nLEGEND_POSITION_Y\t60\nLEGEND_HORIZONTAL\t0\n")
    f.write("LEGEND_SHAPES\t" + "\t".join(["1"] * (len(shown) + 1)) + "\nLEGEND_COLORS\t" + "\t".join([col[g] for g in shown] + ["#D9D9D9"]) + "\nLEGEND_LABELS\t" + "\t".join(shown + ["other genera"]) + "\nDATA\n")
    for n in tips: g = gen.get(n, ""); f.write("%s\t%s\t%s\n" % (n, col.get(g, "#D9D9D9"), g))
G = pd.read_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", dtype=str)
size = G.groupby("clone_group").accession.transform("size"); cs = dict(zip("usr_" + G.accession, size))
with open(f"{out}/dataset_clone_size.txt", "w") as f:
    f.write("DATASET_SIMPLEBAR\nSEPARATOR TAB\nDATASET_LABEL\tGenomes in the clone group (ANI >= 99.9%)\nCOLOR\t#8FA1B3\n")
    f.write("LEGEND_TITLE\tClone group size\nLEGEND_SCALE\t1\nLEGEND_POSITION_X\t60\nLEGEND_POSITION_Y\t400\nLEGEND_HORIZONTAL\t0\nLEGEND_SHAPES\t1\nLEGEND_COLORS\t#8FA1B3\nLEGEND_LABELS\tgenomes in the clone group (ANI >= 99.9%)\nDATA\n")
    for n in tips: f.write("%s\t%d\n" % (n, int(cs.get(n, 1))))
with open(f"{out}/labels.txt", "w") as f:
    f.write("LABELS\nSEPARATOR TAB\nDATA\n")
    for n in tips: f.write("%s\t%s | %s\n" % (n, sp.get(n, "") or gen.get(n, ""), n[4:]))
print("wrote", sorted(os.listdir(out)))
