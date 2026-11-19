import pandas as pd, numpy as np, glob, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
base = "/PATH/TO/virgibacillus_v2"; out = f"{base}/13_figures"; os.makedirs(out, exist_ok=True)
SET = os.environ.get("SET", "main"); NP = int(os.environ.get("NPERM", "300"))
res = sorted(glob.glob(f"{base}/10_pangenome/proteins/OrthoFinder/Results_*"))[-1]
gc = pd.read_csv(f"{res}/Orthogroups/Orthogroups.GeneCount.tsv", sep="\t", index_col=0).drop(columns=["Total"])
un = pd.read_csv(f"{res}/Orthogroups/Orthogroups_UnassignedGenes.tsv", sep="\t", index_col=0, dtype=str)
reps = pd.read_csv(f"{base}/10_pangenome/{'b_reps.tsv' if SET == 'main' else 'b_reps_995.tsv'}", sep="\t", dtype=str)
cols = [a for a in reps.accession if a in gc.columns]; N = len(cols); print("set:", SET, "| genomes:", N)
def lens(a):
    d, cur = {}, None
    for line in open(f"{base}/10_pangenome/proteins/{a}.faa"):
        if line.startswith(">"): cur = line[1:].split()[0]; d[cur] = 0
        else: d[cur] += len(line.strip())
    return d
A0 = gc[cols] > 0; A0 = A0.loc[A0.sum(axis=1) > 0].values.astype(bool)
sg, sl = [], []
for gi, c in enumerate(cols):
    L = lens(c)
    for cell in un[c].dropna():
        for g in str(cell).split(", "):
            g = g.strip()
            if g: sg.append(gi); sl.append(L.get(g, np.nan))
sg = np.array(sg); sl = np.nan_to_num(np.array(sl, float), nan=0.0)
def with_single(mask):
    g = sg[mask]; S = np.zeros((len(g), N), bool); S[np.arange(len(g)), g] = True; return np.vstack([A0, S])
schemes = {"shared orthogroups only": A0, "+ genome-specific genes (all)": with_single(np.ones(len(sg), bool)), "+ genome-specific genes >= 100 aa": with_single(sl >= 100)}
def curves(A, seed=3):
    rng = np.random.default_rng(seed); pan = np.zeros((NP, N)); core = np.zeros((NP, N)); new = np.zeros((NP, N - 1))
    for p in range(NP):
        seen = np.zeros(A.shape[0], bool); inall = np.ones(A.shape[0], bool)
        for k, g in enumerate(rng.permutation(N)):
            col = A[:, g]
            if k > 0: new[p, k - 1] = (col & ~seen).sum()
            seen |= col; inall &= col; pan[p, k] = seen.sum(); core[p, k] = inall.sum()
    return pan, core, new
PAL = ["#7A2F5E", "#D9697A", "#F4B183"]; xs = np.arange(1, N + 1); Ns = np.arange(2, N + 1)
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8), gridspec_kw={"width_ratios": [1.3, 1, 1.1]}); tab = []; core_done = False
for (name, A), col in zip(schemes.items(), PAL):
    pan, core, new = curves(A); lo, hi = np.percentile(pan, [2.5, 97.5], axis=0)
    ax[0].fill_between(xs, lo, hi, color=col, alpha=0.25, lw=0); ax[0].plot(xs, pan.mean(axis=0), color=col, lw=2, label=name)
    if not core_done:
        clo, chi = np.percentile(core, [2.5, 97.5], axis=0); ax[0].fill_between(xs, clo, chi, color="#8FA1B3", alpha=0.35, lw=0)
        ax[0].plot(xs, core.mean(axis=0), color="#3B4252", lw=2, label="core (shared by all genomes sampled)"); core_done = True
    nm = new.mean(axis=0); al = np.array([-np.polyfit(np.log(Ns), np.log(np.maximum(new[p], 1)), 1)[0] for p in range(NP)]); a_lo, a_hi = np.percentile(al, [2.5, 97.5])
    sl_, ic_ = np.polyfit(np.log(Ns), np.log(np.maximum(nm, 1)), 1)
    ax[1].plot(Ns, np.maximum(nm, 1), color=col, lw=1.8, label="%s: alpha = %.2f (%.2f-%.2f)" % (name, al.mean(), a_lo, a_hi)); ax[1].plot(Ns, np.exp(ic_) * Ns ** sl_, color=col, lw=1, ls="--")
    tab.append({"set": SET, "genomes": N, "scheme": name, "families_at_N": round(pan[:, -1].mean()), "alpha": round(al.mean(), 2), "alpha_lo": round(a_lo, 2), "alpha_hi": round(a_hi, 2), "new_in_last_genome": round(new[:, -1].mean(), 1)})
ax[0].set_xlabel("Genomes sampled (random orders, %d permutations)" % NP); ax[0].set_ylabel("Gene families"); ax[0].set_ylim(bottom=0); ax[0].legend(fontsize=7, frameon=False, loc="lower right"); ax[0].set_title("A. Pangenome and core accumulation", fontsize=10)
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("Genomes (N)"); ax[1].set_ylabel("New families added by the N-th genome"); tk = [v for v in (2, 3, 5, 10, 20) if v < N] + [N]; ax[1].set_xticks(tk); ax[1].set_xticklabels([str(v) for v in tk]); ax[1].minorticks_off(); ax[1].legend(fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17)); ax[1].set_title("B. Heaps' law (dashed: fit to the mean)", fontsize=10)
freq = A0.sum(axis=1); c95 = int(np.ceil(0.95 * N)); c15 = int(np.ceil(0.15 * N))
parts = [("core\n(>= 95% of genomes)", int((freq >= c95).sum()), "#3B4252"), ("shell\n(15-95%)", int(((freq >= c15) & (freq < c95)).sum()), "#8C6BB1"), ("cloud\n(< 15%)", int((freq < c15).sum()), "#C07BA0"),
         ("genome-specific\ngenes >= 100 aa", int((sl >= 100).sum()), "#F4B183"), ("genome-specific\ngenes < 100 aa", int((sl < 100).sum()), "#F7C99A")]
ypos = [0, 1, 2, 3.7, 4.7]; tot_og = sum(p[1] for p in parts[:3]); ax[2].barh(ypos, [p[1] for p in parts], color=[p[2] for p in parts], edgecolor="#777777", linewidth=0.4, height=0.75)
for y_, p in zip(ypos, parts): ax[2].text(p[1], y_, "  %d" % p[1] + ("  (%.0f%% of orthogroups)" % (100 * p[1] / tot_og) if y_ < 3 else ""), va="center", fontsize=7.5)
ax[2].set_yticks(ypos); ax[2].set_yticklabels([p[0] for p in parts], fontsize=8); ax[2].invert_yaxis(); ax[2].set_xlim(0, max(p[1] for p in parts) * 1.55); ax[2].set_xlabel("Gene families"); ax[2].set_title("C. Composition (%d genomes)" % N, fontsize=10)
for a in ax: a.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); plt.savefig(f"{out}/fig_pangenome_{SET}.png", dpi=250, bbox_inches="tight"); plt.savefig(f"{out}/fig_pangenome_{SET}.pdf", bbox_inches="tight")
pd.DataFrame(tab).to_csv(f"{out}/pangenome_summary_{SET}.tsv", sep="\t", index=False); print(pd.DataFrame(tab).to_string(index=False)); print([(p[0].replace(chr(10), " "), p[1]) for p in parts]); print("wrote", f"{out}/fig_pangenome_{SET}.png")
