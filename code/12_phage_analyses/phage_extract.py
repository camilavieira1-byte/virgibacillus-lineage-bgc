import os, glob, pandas as pd
base = "/PATH/TO/virgibacillus_v2"; gd = f"{base}/11_genomad"; out = f"{base}/12_phages"
os.makedirs(f"{out}/seqs", exist_ok=True)
rows = []
for d in sorted(glob.glob(f"{gd}/usr_*")):
    name = os.path.basename(d); acc = name[4:]
    fa = f"{d}/{name}_summary/{name}_virus.fna"
    if not os.path.exists(fa): continue
    for rec in open(fa).read().split(">")[1:]:
        head, *seq = rec.split("\n"); sid = head.split()[0]
        if "|provirus_" not in sid: continue
        ctg, rest = sid.split("|provirus_"); s, e = (int(x) for x in rest.split("_")[:2])
        pid = f"{acc}__{ctg}_provirus_{s}_{e}"; sq = "".join(seq)
        open(f"{out}/seqs/{pid}.fna", "w").write(f">{pid}\n{sq}\n")
        rows.append({"pid": pid, "accession": acc, "contig": ctg, "start": s, "end": e, "length": len(sq)})
P = pd.DataFrame(rows)
m = pd.read_csv(f"{gd}/mge_regions.tsv", sep="\t"); m = m[m.kind == "prophage"][["accession", "contig", "start", "end", "in_seq_edge", "virus_score", "n_hallmarks"]]
P = P.merge(m, on=["accession", "contig", "start", "end"], how="left")
P.to_csv(f"{out}/prophages.tsv", sep="\t", index=False)
open(f"{out}/list.txt", "w").write("\n".join(f"{out}/seqs/{p}.fna" for p in P.pid) + "\n")
print("Prophage sequences written:", len(P), "(expected 256) | genomes:", P.accession.nunique(), "| length kb: median %.0f, range %.0f-%.0f" % (P.length.median() / 1e3, P.length.min() / 1e3, P.length.max() / 1e3))
print("rows with an edge flag:", int(P.in_seq_edge.notna().sum()), "| touching a contig end:", int((P.in_seq_edge == 1).sum()))
