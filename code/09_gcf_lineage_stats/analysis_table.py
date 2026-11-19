import pandas as pd, re
base = "/PATH/TO/virgibacillus_v2"
t = pd.read_csv(f"{base}/00_metadata/metadata_curation_template.tsv", sep="\t", dtype=str).fillna("")
s = pd.read_csv(f"{base}/04_gtdbtk_b/gtdb_classification.tsv", sep="\t")
g = pd.read_csv(f"{base}/06_antismash/bgc_per_genome.tsv", sep="\t")
df = g[["accession", "Contig_N50", "n_bgc", "n_edge"]].merge(
    s[["accession", "genus", "species", "Completeness", "Contamination"]], on="accession", how="left").merge(
    t[["accession", "strain", "country", "isolation_source", "habitat_proposed"]], on="accession", how="left")
for c in ["genus", "species", "strain", "country", "isolation_source", "habitat_proposed"]:
    df[c] = df[c].fillna("")

OUT = {"Sediminibacillus", "Terribacillus", "Metabacillus", "Salinibacillus", "Paenibacillus_AQ"}
df["in_A"] = ~df.genus.isin(OUT)
df["in_B"] = df.genus == "Virgibacillus"
df["n50_100k"] = df.Contig_N50 >= 100e3
df["habitat"] = df.habitat_proposed.replace("", "unresolved")

R = {}
def add(region, names):
    for n in names: R[n.lower()] = region
add("Americas", ["USA", "United States", "Puerto Rico", "Brazil", "Mexico", "Canada", "Peru", "Chile", "Argentina", "Colombia"])
add("Europe", ["Spain", "United Kingdom", "UK", "Germany", "France", "Italy", "Greece", "Bulgaria", "Norway", "Netherlands", "Belgium", "Portugal", "Poland", "Russia", "Denmark", "Sweden", "Switzerland", "Austria", "Hungary", "Czech Republic", "Romania", "Ireland"])
add("Africa", ["Algeria", "Tunisia", "Senegal", "Egypt", "Morocco", "South Africa", "Kenya", "Nigeria", "Ethiopia"])
add("Middle East", ["Saudi Arabia", "Turkey", "Iran", "Israel", "Iraq", "Kuwait", "Oman", "United Arab Emirates", "Qatar", "Jordan"])
add("South Asia", ["India", "Pakistan", "Bangladesh", "Sri Lanka", "Nepal"])
add("East Asia", ["China", "Japan", "South Korea", "Republic of Korea", "Korea", "Taiwan", "Mongolia", "North Korea"])
add("Southeast Asia & Oceania", ["Thailand", "Viet Nam", "Vietnam", "Indonesia", "Malaysia", "Philippines", "Singapore", "Micronesia", "Australia", "New Zealand", "Papua New Guinea"])
add("Ocean/Polar", ["Indian Ocean", "Antarctica", "Atlantic Ocean", "Pacific Ocean", "Arctic Ocean", "Southern Ocean", "Arctic", "Atlantic", "Pacific"])
df["region"] = df.country.str.strip().str.lower().map(R).fillna("unresolved")
BAD = {"", "na", "missing", "not applicable", "not collected", "not provided", "unknown", "none"}
un = df[(df.region == "unresolved") & (~df.country.str.strip().str.lower().isin(BAD))].country.value_counts()
print("Countries not placed in a region (tell me where they go):", un.to_dict())

ed = pd.read_csv(f"{base}/05_ani/skani_edges.tsv", sep="\t")
acc = lambda x: re.search(r"(GCA_\d+\.\d)", x).group(1)
ed["a"], ed["b"] = ed.Ref_file.map(acc), ed.Query_file.map(acc)
ed["af"] = ed[["Align_fraction_ref", "Align_fraction_query"]].min(axis=1)
par = {a: a for a in df.accession}
def f(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for r in ed[(ed.ANI >= 99.9) & (ed.af >= 50)].itertuples():
    if r.a in par and r.b in par: par[f(r.a)] = f(r.b)
df["clone_group"] = df.accession.map(f)

rep = df.sort_values(["Completeness", "Contig_N50"], ascending=False).groupby("clone_group").head(1).copy()
for col in ("habitat", "region"):
    mode = df[df[col] != "unresolved"].groupby("clone_group")[col].agg(lambda x: x.value_counts().index[0])
    nmix = df[df[col] != "unresolved"].groupby("clone_group")[col].nunique()
    print("Clone groups with >1 %s label: %d" % (col, int((nmix > 1).sum())))
    rep[col] = rep.clone_group.map(mode).fillna("unresolved")
df.to_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", index=False)
rep.to_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", index=False)

def show(label, d):
    print("\n--- %s: %d clone groups ---" % (label, len(d)))
    print("habitat:", d.habitat.value_counts().to_dict())
    print("region :", d.region.value_counts().to_dict())
show("clade A' (all N50)", rep[rep.in_A])
show("clade A', N50 >= 100 kb", rep[rep.in_A & rep.n50_100k])
show("GTDB Virgibacillus only (B)", rep[rep.in_B])
a = rep[rep.in_A & (rep.habitat != "unresolved") & (rep.region != "unresolved")]
print("\nHabitat x region, clade A', both resolved (%d clone groups):" % len(a))
print(pd.crosstab(a.habitat, a.region).to_string())
