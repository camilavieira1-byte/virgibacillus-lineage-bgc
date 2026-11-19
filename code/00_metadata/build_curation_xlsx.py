import pandas as pd, numpy as np, re, os
from openpyxl import load_workbook
from openpyxl.styles import PatternFill, Font, Alignment
from openpyxl.worksheet.datavalidation import DataValidation
from openpyxl.utils import get_column_letter

base = "/PATH/TO/virgibacillus_v2"; md = f"{base}/00_metadata"

# ---- template + BioSample (NCBI) ----
t = pd.read_csv(f"{md}/metadata_curation_template.tsv", sep="\t", dtype=str)
t = t.drop(columns=[c for c in ["title", "comment"] if c in t.columns])
e = pd.read_csv(f"{md}/biosample_efetch.tsv", sep="\t", dtype=str)
ren = {"title": "bs_title", "comment": "bs_comment", "type-material": "type_material",
       "sample_type": "bs_sample_type", "depth": "bs_depth", "temp": "bs_temp"}
e = e[[c for c in ["accession"] + list(ren) if c in e.columns]].rename(columns=ren)
df = t.merge(e, on="accession", how="left")

# ---- EBI (optional) ----
p_ebi = f"{md}/biosample_ebi.tsv"
if os.path.exists(p_ebi):
    x = pd.read_csv(p_ebi, sep="\t", dtype=str)
    x = x.drop(columns=[c for c in ["biosample", "error"] if c in x.columns])
    x.columns = ["accession" if c == "accession" else "ebi_" + re.sub(r"\W+", "_", c) for c in x.columns]
    df = df.merge(x, on="accession", how="left")

# ---- BioProject / papers ----
b = pd.read_csv(f"{md}/bioproject_papers.tsv", sep="\t", dtype=str).fillna("")
b["accession"] = b.accessions.str.split(";")
b = b.explode("accession").rename(columns={"title": "project_title", "n_genomes": "n_in_project"})
b = b.reindex(columns=["accession", "bioproject", "n_in_project", "project_title", "pmids", "paper_titles", "dois"]).fillna("")
df = df.merge(b, on="accession", how="left")

# ---- clone groups (skani, ANI >= 99.9, AF >= 50) ----
ed = pd.read_csv(f"{base}/05_ani/skani_edges.tsv", sep="\t")
acc = lambda s: re.search(r"(GCA_\d+\.\d)", s).group(1)
ed["a"], ed["b"] = ed.Ref_file.map(acc), ed.Query_file.map(acc)
ed["af"] = ed[["Align_fraction_ref", "Align_fraction_query"]].min(axis=1)
gen = [acc(l) for l in open(f"{base}/05_ani/genome_list.txt")]
par = {g: g for g in gen}
def f(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for r in ed[(ed.ANI >= 99.9) & (ed.af >= 50)].itertuples(): par[f(r.a)] = f(r.b)
groups = {}
for g in gen: groups.setdefault(f(g), []).append(g)
order = sorted(groups.values(), key=lambda c: (-len(c), c[0]))
cg = {g: (f"CG{i+1:03d}", len(c)) for i, c in enumerate(order) for g in c}
df["clone_group"] = df.accession.map(lambda a: cg.get(a, ("not_in_ANI_set", ""))[0])
df["clone_group_size"] = df.accession.map(lambda a: cg.get(a, ("", ""))[1])

# ---- links + manual columns ----
df["biosample_url"] = "https://www.ncbi.nlm.nih.gov/biosample/" + df.biosample.fillna("")
df["bioproject_url"] = np.where(df.bioproject.fillna("") != "", "https://www.ncbi.nlm.nih.gov/bioproject/" + df.bioproject.fillna(""), "")
df["doi_url"] = df.dois.fillna("").map(lambda s: ("https://doi.org/" + s.split(";")[0]) if s else "")
df = df.rename(columns={"lat": "lat_biosample", "lon": "lon_biosample"})
MANUAL = ["lat_manual", "lon_manual", "location_text_manual", "coords_precision", "source_ref", "habitat_final", "notes"]
for c in MANUAL: df[c] = ""

df["n_in_project"] = pd.to_numeric(df.n_in_project, errors="coerce").fillna(0)
df = df.sort_values(["n_in_project", "bioproject", "clone_group", "accession"], ascending=[False, True, True, True])

front = ["accession", "strain", "type_material", "organism", "pass_qc", "clone_group", "clone_group_size", "bioproject",
         "country", "geo_loc_name", "isolation_source", "bs_title", "lat_biosample", "lon_biosample"]
rest = [c for c in df.columns if c not in front + MANUAL]
df = df[front + MANUAL + rest]
df.to_csv(f"{md}/metadata_curation.tsv", sep="\t", index=False)

# ---- Excel ----
out = f"{md}/metadata_curation.xlsx"
with pd.ExcelWriter(out, engine="openpyxl") as w:
    df.to_excel(w, sheet_name="curation", index=False)
    b.drop_duplicates("bioproject").to_excel(w, sheet_name="projects", index=False)
    pd.DataFrame({"README": [
        "Fill ONLY the yellow columns in sheet 'curation'. Everything else is context.",
        "lat_manual / lon_manual: decimal degrees, WGS84, dot as decimal separator (-23.5, 46.6). South and West are negative.",
        "coords_precision: exact | locality_approx | region_approx | country_only | none",
        "source_ref: DOI, PMID, BioSample/BioProject URL, culture collection page, or 'species description: <ref>'.",
        "habitat_final: pick from the dropdown. Leave blank only if truly unknown (then choose 'unknown').",
        "Rows in the same clone_group (CG###) are near-identical genomes: curate one and copy to the others if they share a source.",
        "Rows with lat_biosample already filled: check them, and fill coords_precision. Do not overwrite unless wrong (explain in notes).",
        "Type strains (type_material): coordinates usually come from the original species description, not the BioSample.",
        "Do not insert or delete columns, and do not rename the headers (the ingest script reads them by name).",
    ]}).to_excel(w, sheet_name="README", index=False)

wb = load_workbook(out); ws = wb["curation"]
hdr = [c.value for c in ws[1]]
yellow = PatternFill("solid", fgColor="FFF2CC"); head = PatternFill("solid", fgColor="D9E1F2"); headm = PatternFill("solid", fgColor="FFD966")
for j, name in enumerate(hdr, 1):
    L = get_column_letter(j)
    ws[f"{L}1"].font = Font(bold=True); ws[f"{L}1"].fill = headm if name in MANUAL else head
    vals = [len(str(ws.cell(r, j).value or "")) for r in range(1, min(ws.max_row, 80) + 1)]
    ws.column_dimensions[L].width = max(10, min(max(vals) + 2, 45))
    if name in MANUAL:
        for r in range(2, ws.max_row + 1): ws.cell(r, j).fill = yellow
    if name.endswith("_url"):
        for r in range(2, ws.max_row + 1):
            c = ws.cell(r, j)
            if c.value: c.hyperlink = str(c.value); c.font = Font(color="0563C1", underline="single")
ws.freeze_panes = "C2"; ws.auto_filter.ref = ws.dimensions
col = lambda n: get_column_letter(hdr.index(n) + 1)
dv1 = DataValidation(type="list", formula1='"constructed,soil,saline_lake,saltern,marine,mangrove,host_associated,food,other,unknown"', allow_blank=True)
dv2 = DataValidation(type="list", formula1='"exact,locality_approx,region_approx,country_only,none"', allow_blank=True)
ws.add_data_validation(dv1); ws.add_data_validation(dv2)
dv1.add(f"{col('habitat_final')}2:{col('habitat_final')}{ws.max_row}")
dv2.add(f"{col('coords_precision')}2:{col('coords_precision')}{ws.max_row}")
ws2 = wb["README"]; ws2.column_dimensions["A"].width = 140
wb.save(out)

print("Rows:", len(df), "| pass QC:", int((df.pass_qc == "True").sum()))
print("Clone groups with >1 genome:", sum(1 for c in groups.values() if len(c) > 1), "| distinct clone groups:", len(groups))
print("Rows already with lat/lon (BioSample):", int(df.lat_biosample.notna().sum()))
print("Rows needing coordinates:", int(df.lat_biosample.isna().sum()))
print("Wrote:", out)
