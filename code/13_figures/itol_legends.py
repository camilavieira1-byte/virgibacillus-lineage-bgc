import pandas as pd, re, os
base = "/PATH/TO/virgibacillus_v2"; d = f"{base}/13_figures/itol"
def patch(fn, color=None, extra=None):
    L = open(f"{d}/{fn}").read().split("\n"); out = []
    for l in L:
        if color and l.startswith("COLOR\t"): l = "COLOR\t" + color
        if l == "DATA" and extra: out += extra
        out.append(l)
    open(f"{d}/{fn}", "w").write("\n".join(out))
classes = ["NRPS", "NRP-metallophore", "RiPP-like", "lanthipeptide-class-ii", "betalactone"]
colors = ["#7A2F5E", "#D9697A", "#F4B183", "#8C6BB1", "#4C78A8"]
def legend(title, shapes, cols, labels, y):
    return ["LEGEND_TITLE\t" + title, "LEGEND_SCALE\t1", "LEGEND_POSITION_X\t60", "LEGEND_POSITION_Y\t%d" % y, "LEGEND_HORIZONTAL\t0",
            "LEGEND_SHAPES\t" + "\t".join(shapes), "LEGEND_COLORS\t" + "\t".join(cols), "LEGEND_LABELS\t" + "\t".join(labels)]
patch("dataset_bgc_classes.txt", extra=legend("BGC classes (circle = class present in the genome)", ["2"] * 5, colors, classes, 160))
PRO, BGC = "#3B4252", "#8FA1B3"
patch("dataset_prophages.txt", color=PRO, extra=["DATASET_SCALE\t2\t4\t6\t8"] + legend("Prophages per genome (geNomad)", ["1"], [PRO], ["bar length = number of prophage regions (scale lines at 2, 4, 6, 8)"], 340))
patch("dataset_bgc_regions.txt", color=BGC, extra=["DATASET_SCALE\t4\t8\t12\t16"] + legend("BGC regions per genome (antiSMASH 8.0.4)", ["1"], [BGC], ["bar length = number of regions, terpene-precursor excluded (scale lines at 4, 8, 12, 16)"], 420))
G = pd.read_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", dtype=str).fillna("").set_index("accession")
n50 = pd.to_numeric(G.Contig_N50, errors="coerce")
tips = [l.split("\t")[0] for l in open(f"{d}/labels.txt").read().split("\n")[3:] if l.strip()]
def band(x):
    if x != x: return "#FFFFFF", "unknown"
    return ("#252525", "N50 < 100 kb") if x < 100e3 else (("#969696", "N50 100-500 kb") if x < 500e3 else ("#D9D9D9", "N50 >= 500 kb"))
with open(f"{d}/dataset_n50.txt", "w") as f:
    f.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tAssembly N50\nCOLOR\t#969696\n")
    f.write("\n".join(legend("Assembly contiguity (N50)", ["1", "1", "1"], ["#252525", "#969696", "#D9D9D9"], ["N50 < 100 kb", "N50 100-500 kb", "N50 >= 500 kb"], 520)) + "\nDATA\n")
    for t in tips:
        a = re.sub(r"^usr_", "", t); c, l = band(n50.get(a, float("nan"))); f.write(f"{t}\t{c}\t{l}\n")
print(sorted(os.listdir(d)))
