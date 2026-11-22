import pandas as pd, numpy as np, os, collections
from Bio import Phylo
from scipy.optimize import minimize
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
base = "/PATH/TO/virgibacillus_v2"; out = f"{base}/13_figures"; os.makedirs(out, exist_ok=True)
MODEL = os.environ.get("MODEL", "ARD"); CLASSES = os.environ.get("CLASSES", ""); TIPS_ONLY = os.environ.get("TIPS_ONLY") == "1"
MIN_N50 = float(os.environ.get("MIN_N50", "0")); TAG = os.environ.get("TAG", "all")
cg = pd.read_csv(f"{base}/09_matrix/{os.environ.get('TABLE','analysis_table_clonegroups.tsv')}", sep="\t", dtype=str).fillna("")
cg = cg[cg.in_A == "True"].copy(); cg["n50"] = pd.to_numeric(cg.Contig_N50, errors="coerce"); cg = cg[cg.n50 >= MIN_N50]
if os.environ.get("ONLY_B") == "1": cg = cg[cg.in_B == "True"]
gt = pd.read_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", dtype=str).fillna(""); info = cg.set_index("accession")
t = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick"); tipd = {x.name: x for x in t.get_terminals()}
og = ["usr_" + a for a in gt[gt.genus.str.startswith("Paenibacillus")].accession if "usr_" + a in tipd]
if og: t.root_with_outgroup(*[tipd[nm] for nm in og])
else: print("WARNING: no outgroup tip found; tree not rooted")
keep = {"usr_" + a for a in cg.accession}
for x in list(t.get_terminals()):
    if x.name not in keep: t.prune(x)
while len(t.root.clades) == 1: t.root = t.root.clades[0]
tips = t.get_terminals(); n = len(tips); print("tips in tree: %d (expected %d) | N50 filter: >= %.0f bp" % (n, len(keep), MIN_N50))
post = list(t.find_clades(order="postorder")); pre = list(t.find_clades(order="preorder"))
yv = {id(x): i for i, x in enumerate(tips)}; h = {}
for c in post: h[id(c)] = 0 if c.is_terminal() else 1 + max(h[id(k)] for k in c.clades)
H = h[id(t.root)]; xv = {id(c): H - h[id(c)] for c in post}
for c in post:
    if not c.is_terminal(): yv[id(c)] = float(np.mean([yv[id(k)] for k in c.clades]))
bl = {id(c): max(c.branch_length or 0.0, 1e-6) for c in post}
for c in post:
    c.branch_length = None if c is t.root else bl[id(c)]; c.confidence = None
    if not c.is_terminal(): c.name = None
Phylo.write(t, f"{out}/asr_tree_{TAG}.nwk", "newick")
Bset = {"usr_" + a for a in cg[cg.in_B == "True"].accession}
grp = lambda nm: info.species.get(nm[4:], "") or info.genus.get(nm[4:], "")
label = lambda nm: grp(nm) + "  " + nm[4:]
reg = pd.read_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t"); tk = reg.assign(tok=reg["product"].str.split("+")).explode("tok")
MIN_GENOMES = int(os.environ.get("MIN_GENOMES", "8"))
_cnt = tk[tk.accession.isin(set(cg.accession))].groupby("tok").accession.nunique()
classes = sorted(_cnt[_cnt >= MIN_GENOMES].index.tolist())
print("classes present in >= %d of the %d representatives: %d" % (MIN_GENOMES, len(set(cg.accession)), len(classes)))
for _c in classes: print("   %-26s %d genomes" % (_c, _cnt[_c]))
_pal = ["#7A2F5E", "#D9697A", "#F4B183", "#8C6BB1", "#4C78A8", "#E8999B", "#8FA1B3",
        "#2F6E5E", "#B5651D", "#5A3E8C", "#7A1F2F", "#3C7A2F", "#C2A24C", "#4F6D7A"]
colors = [_pal[i % len(_pal)] for i in range(len(classes))]
run_classes = [c for c in classes if (not CLASSES) or c in CLASSES.split(",")]
pres = {c: set(tk[tk.tok == c].accession) for c in classes}
pd.DataFrame({c: {x.name: int(x.name[4:] in pres[c]) for x in tips} for c in classes}).to_csv(f"{out}/asr_traits_{TAG}.tsv", sep="\t")
open(f"{out}/asr_Btips_{TAG}.txt", "w").write("\n".join(sorted(Bset)))

def transition(a, b, tl):
    s = a + b; e = np.exp(-s * tl)
    return np.array([[(b + a * e) / s, a * (1 - e) / s], [b * (1 - e) / s, (a + b * e) / s]])
def downpass(a, b, st):
    Pm = {id(c): transition(a, b, bl[id(c)]) for c in post if c is not t.root}; down = {}
    for c in post:
        if c.is_terminal(): down[id(c)] = np.array([1.0, 0.0]) if st[c.name] == 0 else np.array([0.0, 1.0])
        else:
            v = np.ones(2)
            for k in c.clades: v = v * (Pm[id(k)] @ down[id(k)])
            down[id(c)] = v
    return Pm, down
def nll(p, st):
    a, b = np.exp(np.clip(p, -9, 9)); Pm, down = downpass(a, b, st)
    return -np.log((0.5 * down[id(t.root)]).sum() + 1e-300)
def fit(st):
    best = None
    for x0 in ([0, 0], [-2, 1], [1, -2], [3, 3]):
        r = minimize(nll, x0, args=(st,), method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 800})
        if best is None or r.fun < best.fun: best = r
    return np.exp(np.clip(best.x, -9, 9))
def fit_er(st):
    best = None
    for x0 in ([0.0], [3.0], [-2.0]):
        r = minimize(lambda p, s_: nll([p[0], p[0]], s_), x0, args=(st,), method="Nelder-Mead")
        if best is None or r.fun < best.fun: best = r
    q = float(np.exp(np.clip(best.x[0], -9, 9))); return np.array([q, q])
def marginal(a, b, st):
    Pm, down = downpass(a, b, st); up = {id(t.root): np.array([0.5, 0.5])}
    for c in pre:
        for k in c.clades:
            ex = up[id(c)].copy()
            for s2 in c.clades:
                if s2 is not k: ex = ex * (Pm[id(s2)] @ down[id(s2)])
            up[id(k)] = ex @ Pm[id(k)]
    res = {}
    for c in post:
        p = down[id(c)] * up[id(c)]; res[id(c)] = p / p.sum()
    return res
def fitch(st):
    sets, ch = {}, 0
    for c in post:
        if c.is_terminal(): sets[id(c)] = {st[c.name]}
        else:
            kids = [sets[id(k)] for k in c.clades]; inter = set.intersection(*kids)
            if inter: sets[id(c)] = inter
            else: sets[id(c)] = set.union(*kids); ch += 1
    return ch

bt = [x for x in tips if x.name in Bset]; mB = t.common_ancestor(bt) if len(bt) > 1 else None
spB = collections.defaultdict(list)
for x in bt: spB[grp(x.name)].append(x)
par = {id(k): c for c in pre for k in c.clades}; results, rows, at_bound = {}, [], {}
for cl in run_classes:
    st = {x.name: int(x.name[4:] in pres[cl]) for x in tips}; k = sum(st.values())
    if k in (0, n): print("\n==", cl, ": constant across tips, skipped"); continue
    a, b = fit_er(st) if MODEL == "ER" else fit(st); pp = marginal(a, b, st); results[cl] = (st, pp); gains, losses = [], []
    at_bound[cl] = bool(max(a, b) >= 0.99 * np.exp(9) or min(a, b) <= 1.01 * np.exp(-9))
    for c in pre:
        if c is t.root: continue
        p0, p1 = pp[id(par[id(c)])][1], pp[id(c)][1]
        if p0 < 0.5 <= p1: gains.append(c)
        if p0 >= 0.5 > p1: losses.append(c)
    pB = pp[id(mB)][1] if mB is not None else np.nan
    print("\n== %s: present in %d of %d tips | rates per unit branch length: gain %.2f, loss %.2f | P(present) at root %.2f, at the Virgibacillus ancestor %.2f" % (cl, k, n, a, b, pp[id(t.root)][1], pB))
    print("   inferred gains: %d | inferred losses: %d | minimum changes (Fitch parsimony): %d" % (len(gains), len(losses), fitch(st)))
    for c in sorted(gains, key=lambda c: -len(c.get_terminals()))[:8]:
        comp = collections.Counter(grp(x.name) for x in c.get_terminals()).most_common(3)
        print("   gain on a branch leading to %d tips (P %.2f -> %.2f): %s" % (len(c.get_terminals()), pp[id(par[id(c)])][1], pp[id(c)][1], "; ".join("%s x%d" % z for z in comp)))
    for sp, xs in spB.items():
        if len(xs) > 1: print("   species ancestor %-28s P(present) = %.2f" % (sp, pp[id(t.common_ancestor(xs))][1]))
    rows.append({"class": cl, "tips_present": k, "tips": n, "gain_rate": a, "loss_rate": b, "P_root": pp[id(t.root)][1], "P_Virgibacillus_ancestor": pB, "gains": len(gains), "losses": len(losses), "fitch_changes": fitch(st), "rates_at_bound": at_bound[cl]})
pd.DataFrame(rows).to_csv(f"{out}/ancestral_summary_{TAG}.tsv", sep="\t", index=False)

cl_all = list(results.keys()); internal = [c for c in post if not c.is_terminal()]; by = sorted(yv[id(x)] for x in bt)
NCOL = int(os.environ.get("NCOL", "5"))
chunks = [cl_all[i:i + NCOL] for i in range(0, len(cl_all), NCOL)]
for _ci, cl_list in enumerate(chunks, 1):
    sfx = "" if len(chunks) == 1 else "_part%d" % _ci
    fig, axs = plt.subplots(1, len(cl_list), figsize=(3.4 * len(cl_list) + 4.5, max(11, 0.17 * n + 2)), sharey=True, squeeze=False); axs = axs[0]
    for ax, cl in zip(axs, cl_list):
        col = colors[classes.index(cl)]; st, pp = results[cl]; cmap = LinearSegmentedColormap.from_list("c", ["#FFFFFF", col])
        if by: ax.axhspan(min(by) - 0.5, max(by) + 0.5, color="#F6E8EF", zorder=0)
        for c in internal:
            ys = [yv[id(k)] for k in c.clades]; ax.plot([xv[id(c)]] * 2, [min(ys), max(ys)], color="#BBBBBB", lw=0.8, zorder=1)
            for k in c.clades:
                on = (pp[id(k)][1] >= 0.5) and not TIPS_ONLY
                ax.plot([xv[id(c)], xv[id(k)]], [yv[id(k)]] * 2, color=col if on else "#BBBBBB", lw=1.4 if on else 0.8, zorder=2)
        if not TIPS_ONLY: ax.scatter([xv[id(c)] for c in internal], [yv[id(c)] for c in internal], c=[pp[id(c)][1] for c in internal], cmap=cmap, vmin=0, vmax=1, s=14, edgecolors="#444444", linewidths=0.3, zorder=3)
        if mB is not None and not TIPS_ONLY: ax.annotate("P = %.2f" % pp[id(mB)][1], (xv[id(mB)], yv[id(mB)]), xytext=(-5, 8), textcoords="offset points", ha="right", fontsize=7, color="#444444", zorder=4)
        for x in tips:
            on = st[x.name] == 1
            ax.scatter([H + 0.5], [yv[id(x)]], marker="s", s=24, facecolors=col if on else "white", edgecolors=col if on else "#BBBBBB", linewidths=0.6, zorder=3)
        ax.set_xlim(-0.5, H + 1.0); ax.set_ylim(n - 0.5, -0.5); ax.set_xticks([]); ax.set_title("%s\n%d of %d tips%s" % (cl, sum(st.values()), n, " (tip pattern only)" if TIPS_ONLY else (" - rates at bound" if at_bound.get(cl) else "")), fontsize=9)
        for s_ in ax.spines.values(): s_.set_visible(False)
    axs[0].set_yticks(range(n)); axs[0].set_yticklabels([label(x.name) for x in tips], fontsize=6)
    for ax in axs[:-1]: ax.tick_params(axis="y", left=False, labelleft=False)
    axs[-1].tick_params(axis="y", left=False, labelleft=False, right=True, labelright=True, length=0)
    if by: axs[0].text(-0.45, (min(by) + max(by)) / 2, "Virgibacillus (GTDB)", rotation=90, va="center", ha="center", fontsize=8, color="#7A2F5E")
    fig.text(0.01, 0.004, "Cladogram (branch lengths not to scale). Square = class observed in the genome. Node shade = marginal probability of presence (ML, " + ("equal-rates" if MODEL == "ER" else "two-rate") + " Mk model, flat root prior); coloured branches: presence inferred (P >= 0.5).", fontsize=7)
    (fig.texts[-1].set_text("Cladogram (branch lengths not to scale). Square = class observed in the genome. Ancestral states are not shown: rate estimates for these classes were unstable.") if TIPS_ONLY else None); plt.tight_layout(rect=(0, 0.015, 1, 1)); plt.savefig(f"{out}/fig_bgc_ancestral_{TAG}{sfx}.png", dpi=200); plt.savefig(f"{out}/fig_bgc_ancestral_{TAG}{sfx}.pdf")
    plt.close(fig)
    print("wrote", f"{out}/fig_bgc_ancestral_{TAG}{sfx}.png")
