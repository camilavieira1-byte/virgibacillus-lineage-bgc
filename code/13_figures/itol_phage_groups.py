import pandas as pd, re
base = "/PATH/TO/virgibacillus_v2"; d = f"{base}/13_figures/itol"
P = pd.read_csv(f"{base}/12_phages/prophages.tsv", sep="\t")
C = pd.read_csv(f"{base}/12_phages/clusters_phabox_aai.tsv", sep="\t", header=None, names=["pid", "c"])
D = P.merge(C, on="pid")
tips = [l.split("\t")[0] for l in open(f"{d}/labels.txt").read().split("\n")[3:] if l.strip()]
acc = {t: re.sub(r"^usr_", "", t) for t in tips}
sub = D[D.accession.isin(set(acc.values()))]
cnt = sub.groupby("c").accession.nunique().sort_values(ascending=False)
groups = [g for g in cnt.index if cnt[g] >= 3][:3]
cols = ["#2A9D8F", "#E9C46A", "#264653"][:len(groups)]
med = {g: sub[sub.c == g].length.median() / 1e3 for g in groups}
print("groups drawn:", {g: "%d genomes, median %.0f kb" % (cnt[g], med[g]) for g in groups})
labels = ["Prophage AAI group %s (%d genomes, median %.0f kb)" % (g.replace("group_", ""), cnt[g], med[g]) for g in groups]
with open(f"{d}/dataset_phage_groups.txt", "w") as f:
    f.write("DATASET_BINARY\nSEPARATOR TAB\nDATASET_LABEL\tProphage groups\nCOLOR\t#264653\n")
    f.write("FIELD_SHAPES\t" + "\t".join(["1"] * len(groups)) + "\nFIELD_LABELS\t" + "\t".join(labels) + "\nFIELD_COLORS\t" + "\t".join(cols) + "\n")
    f.write("LEGEND_TITLE\tProphage groups (PhaBOX2 AAI, genus-level)\nLEGEND_SCALE\t1\nLEGEND_POSITION_X\t60\nLEGEND_POSITION_Y\t600\nLEGEND_HORIZONTAL\t0\n")
    f.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(groups)) + "\nLEGEND_COLORS\t" + "\t".join(cols) + "\nLEGEND_LABELS\t" + "\t".join(labels) + "\nDATA\n")
    for t in tips:
        have = set(sub[sub.accession == acc[t]].c)
        f.write(t + "\t" + "\t".join("1" if g in have else "-1" for g in groups) + "\n")
print("wrote", f"{d}/dataset_phage_groups.txt")
