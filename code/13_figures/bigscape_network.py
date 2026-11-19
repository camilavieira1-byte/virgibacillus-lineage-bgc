import pandas as pd, numpy as np, glob, os, collections
import networkx as nx
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
base = "/PATH/TO/virgibacillus_v2"; out = f"{base}/13_figures"; CUT = os.environ.get("CUT", "0.5")
run = glob.glob(f"{base}/07_bigscape/out_glocal/output_files/*_c{CUT}")[0]
files = sorted(glob.glob(f"{run}/*/*_c{CUT}.network"))
E = pd.concat([pd.read_csv(f, sep="\t", usecols=["GBK_a", "GBK_b", "distance"]) for f in files], ignore_index=True)
E = E[E.distance <= float(CUT)]
rec = pd.read_csv(f"{base}/09_matrix/gcf_records_c{CUT}.tsv", sep="\t")
st = pd.read_csv(f"{base}/06_antismash/mibig_strict_per_region.tsv", sep="\t")[["gbk", "strict"]].drop_duplicates("gbk")
rec = rec.merge(st, on="gbk", how="left"); rec["strict"] = rec.strict.fillna("no hit")
def grp(p):
    if "NRP-metallophore" in p: return "NRP-metallophore"
    if "NRPS" in p or "T1PKS" in p or "transAT" in p: return "NRPS / T1PKS / transAT-PKS"
    if "T3PKS" in p: return "T3PKS"
    if p == "terpene-precursor": return "terpene-precursor"
    if "terpene" in p: return "terpene"
    if "ectoine" in p: return "ectoine"
    if "betalactone" in p: return "beta-lactone"
    if any(k in p for k in ("RiPP", "lanthipeptide", "lassopeptide", "triceptide", "RRE", "azole", "thiopeptide", "LAP")): return "RiPP"
    return "other"
rec["grp"] = rec["product"].map(grp)
COL = {"T3PKS": "#E8999B", "terpene": "#8FA1B3", "terpene-precursor": "#D9D9D9", "ectoine": "#8C6BB1", "RiPP": "#F4B183",
       "NRPS / T1PKS / transAT-PKS": "#7A2F5E", "NRP-metallophore": "#D9697A", "beta-lactone": "#4C78A8", "other": "#EBD0DE"}
G = nx.Graph(); G.add_nodes_from(rec.gbk)
missing = (set(E.GBK_a) | set(E.GBK_b)) - set(rec.gbk)
G.add_edges_from([(a, b) for a, b in zip(E.GBK_a, E.GBK_b) if a in G and b in G])
comps = sorted(nx.connected_components(G), key=len, reverse=True)
multi = [c for c in comps if len(c) >= 2]; solo = [next(iter(c)) for c in comps if len(c) == 1]
print("cutoff %s | nodes %d | edges %d | components with >= 2 regions: %d | isolated regions: %d (singletons in my table: %d) | edge endpoints missing from my table: %d" % (
    CUT, G.number_of_nodes(), G.number_of_edges(), len(multi), len(solo), int(rec.singleton.sum()), len(missing)))
rad = lambda n: 1.5 + 0.7 * np.sqrt(n)
W = 1.1 * np.sqrt(sum((2 * rad(len(c)) + 0.8) ** 2 for c in multi))
pos, centers = {}, []; x = y = rowh = 0.0
for c in multi:
    n = len(c); r = rad(n); sub = G.subgraph(c)
    if n == 2:
        a, b = list(c); loc = {a: np.array([-0.7, 0.0]), b: np.array([0.7, 0.0])}
    else:
        lay = nx.spring_layout(sub, seed=3, k=1.6 / np.sqrt(n), iterations=250); arr = np.array(list(lay.values())); arr = arr - arr.mean(0)
        m = np.abs(arr).max() or 1.0; loc = {k_: v / m * r for k_, v in zip(lay.keys(), arr)}
    bw = 2 * r + 0.8
    if x + bw > W and x > 0: x = 0.0; y += rowh + 0.6; rowh = 0.0
    cx, cy = x + r + 0.4, y + 1.2 + r
    for k_, v in loc.items(): pos[k_] = (cx + v[0], cy + v[1])
    centers.append((c, cx, cy, r)); x += bw; rowh = max(rowh, 2 * r + 1.2)
y += rowh + 1.5; cols = max(int(W / 1.0), 10); solo_y = y
for i, nd in enumerate(solo): pos[nd] = ((i % cols) * 1.0 + 0.5, y + 0.8 + (i // cols) * 1.0)
xs = np.array([p[0] for p in pos.values()]); ys = np.array([p[1] for p in pos.values()])
dw = xs.max() - xs.min() + 3; dh = ys.max() - ys.min() + 4
main_w, leg_w = 12.0, 3.2; main_h = min(22, max(8, main_w * dh / dw))
ppu = min(main_w * 72 / dw, main_h * 72 / dh) * 0.95
size_of = {}
for c, cx, cy, r in centers:
    d = float(np.clip(1.18 * r / np.sqrt(len(c)), 0.45, 0.9))
    for nd in c: size_of[nd] = (d * ppu) ** 2
for nd in solo: size_of[nd] = (0.55 * ppu) ** 2
grp_of, st_of, fam_of, gen_of = dict(zip(rec.gbk, rec.grp)), dict(zip(rec.gbk, rec.strict)), dict(zip(rec.gbk, rec.family)), dict(zip(rec.gbk, rec.genome))
nodes = list(pos)
ec = ["#000000" if st_of[n_].startswith("known") else ("#4A4A4A" if st_of[n_].startswith("short") else "#BBBBBB") for n_ in nodes]
lw = [1.6 if st_of[n_].startswith("known") else (1.0 if st_of[n_].startswith("short") else 0.4) for n_ in nodes]
fig, (ax, axl) = plt.subplots(1, 2, figsize=(main_w + leg_w, main_h), gridspec_kw={"width_ratios": [main_w, leg_w]})
ax.add_collection(LineCollection([(pos[a], pos[b]) for a, b in G.edges()], colors="#9a9a9a", linewidths=0.4, alpha=0.3, zorder=1))
ax.scatter([pos[n_][0] for n_ in nodes], [pos[n_][1] for n_ in nodes], s=[size_of[n_] for n_ in nodes], c=[COL[grp_of[n_]] for n_ in nodes], edgecolors=ec, linewidths=lw, zorder=3)
for c, cx, cy, r in centers:
    if len(c) >= 12:
        top = collections.Counter(fam_of[n_] for n_ in c).most_common(1)[0][0]
        ax.text(cx, cy - r - 0.3, "%s, n=%d" % (top, len(c)), ha="center", va="bottom", fontsize=7)
ax.text(0.2, solo_y + 0.2, "Unconnected regions (n = %d)" % len(solo), fontsize=8, va="bottom")
ax.set_xlim(xs.min() - 1.5, xs.max() + 1.5); ax.set_ylim(ys.max() + 2, ys.min() - 2); ax.set_aspect("equal"); ax.axis("off"); axl.axis("off")
cnt = rec.grp.value_counts()
h1 = [Line2D([0], [0], marker="o", color="w", markerfacecolor=v, markeredgecolor="#777777", markeredgewidth=0.4, markersize=9, label="%s (%d)" % (k_, cnt.get(k_, 0))) for k_, v in COL.items()]
leg1 = axl.legend(handles=h1, loc="upper left", frameon=False, fontsize=8, title="Product class (fill), regions", title_fontsize=8); axl.add_artist(leg1)
h2 = [Line2D([0], [0], marker="o", color="w", markerfacecolor="#EEEEEE", markeredgecolor=e_, markeredgewidth=w_, markersize=9, label=l_)
      for e_, w_, l_ in (("#000000", 1.6, "known (>= 50% similarity, >= 3 genes)"), ("#4A4A4A", 1.0, "short match (>= 50%, < 3 genes)"), ("#BBBBBB", 0.4, "weak or no MIBiG hit"))]
axl.legend(handles=h2, loc="upper left", bbox_to_anchor=(0.0, max(0.2, 1 - (len(h1) + 2.4) * 0.30 / main_h)), frameon=False, fontsize=8, title="MIBiG status (node outline)", title_fontsize=8)
fig.text(0.01, 0.005, "BiG-SCAPE 2.0.3 network, cutoff %s (glocal alignment). Edges: distance <= %s. Each connected component is laid out separately (force-directed); positions between components are arbitrary." % (CUT, CUT), fontsize=7)
plt.tight_layout(rect=(0, 0.015, 1, 1)); plt.savefig(f"{out}/fig_bigscape_network_c{CUT}.png", dpi=250); plt.savefig(f"{out}/fig_bigscape_network_c{CUT}.pdf")
rows = []
for i, (c, cx, cy, r) in enumerate(centers, 1):
    g_ = collections.Counter(grp_of[n_] for n_ in c); f_ = collections.Counter(fam_of[n_] for n_ in c)
    rows.append({"component": i, "regions": len(c), "genomes": len({gen_of[n_] for n_ in c}), "families": len(f_), "top_family": f_.most_common(1)[0][0],
                 "top_group": g_.most_common(1)[0][0], "pct_known": round(100 * np.mean([st_of[n_].startswith("known") for n_ in c]), 1)})
pd.DataFrame(rows).to_csv(f"{out}/network_components_c{CUT}.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).head(10).to_string(index=False)); print("wrote", f"{out}/fig_bigscape_network_c{CUT}.png")
