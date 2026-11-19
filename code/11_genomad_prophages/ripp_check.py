import pandas as pd, re
from Bio import SeqIO
base = "/PATH/TO/virgibacillus_v2"
reg = pd.read_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t")
mge = pd.read_csv(f"{base}/11_genomad/mge_regions.tsv", sep="\t"); pro = mge[mge.kind == "prophage"]
targets = [("GCA_000445495.1", "AUQA01000008.1", 111605), ("GCA_000723585.1", "CCDP010000003.1", 56707), ("GCA_001895305.1", "CP017762.1", 4594328),
           ("GCA_042643265.1", "JBHRRZ010000015.1", 165436), ("GCA_055329585.1", "JBUZWC010000007.1", 19114)]
def orig_start(p):
    for i, l in enumerate(open(p)):
        m = re.search(r"Orig\. start\s*::\s*(\d+)", l)
        if m: return int(m.group(1)) + 1
        if i > 30: return None
for acc, ctg, bs in targets:
    cand = reg[(reg.accession == acc) & (reg.contig == ctg)]
    for f in cand.file:
        p = f"{base}/06_antismash/{acc}/{f}"
        if orig_start(p) != bs: continue
        r = cand[cand.file == f].iloc[0]; ph = pro[(pro.accession == acc) & (pro.contig == ctg)]
        print("\n==", acc, ctg, bs, "|", r["product"], "| contig_edge:", r.contig_edge, "| prophages on this contig (start, end):", ph[["start", "end"]].values.tolist())
        for rec in SeqIO.parse(p, "genbank"):
            for ft in rec.features:
                if ft.type == "CDS": q = ft.qualifiers; print("  ", q.get("gene_kind", [""])[0].ljust(14), q.get("product", [""])[0][:75])
