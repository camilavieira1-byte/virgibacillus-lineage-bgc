import pandas as pd, re, numpy as np
base = "/PATH/TO/virgibacillus_v2"
m = pd.read_csv(f"{base}/00_metadata/metadata_raw_flagged.tsv", sep="\t")
qc = pd.read_csv(f"{base}/03_qc/checkm2/quality_report.tsv", sep="\t")[["Name","Completeness","Contamination"]]
qc = qc.rename(columns={"Name":"accession"})
m = m[m.is_mag=="isolate"].merge(qc, on="accession", how="inner")
m["pass_qc"] = (m.Completeness>=90) & (m.Contamination<5)

BAD = {"na","missing","not applicable","not collected","not provided","unknown","none","nan",""}
def clean(x):
    x = str(x).strip()
    return "" if x.lower() in BAD else x

def parse_latlon(s):
    s = clean(s)
    if not s: return (np.nan, np.nan)
    mm = re.match(r"^\s*([\d.\-]+)\s*([NS])\s*[, ]\s*([\d.\-]+)\s*([EW])\s*$", s, re.I)
    if mm:
        lat = float(mm.group(1)) * (1 if mm.group(2).upper()=="N" else -1)
        lon = float(mm.group(3)) * (1 if mm.group(4).upper()=="E" else -1)
        return (lat, lon)
    mm = re.match(r"^\s*(-?[\d.]+)\s*[, ]\s*(-?[\d.]+)\s*$", s)
    if mm: return (float(mm.group(1)), float(mm.group(2)))
    return (np.nan, np.nan)

m[["lat","lon"]] = m.lat_lon.apply(lambda s: pd.Series(parse_latlon(s)))
m["country"] = m.geo_loc_name.apply(lambda s: clean(s).split(":")[0].strip())

RULES = [  # order matters: first match wins
 ("constructed",     r"clean ?room|spacecraft|assembly facility|jpl|floor"),
 ("host_associated", r"stool|feces|faec|gut|intestin|human|host|skin|oral"),
 ("food",            r"kimchi|ferment|food|fish sauce|sauce|cheese|soy|miso|salted|jeotgal"),
 ("mangrove",        r"mangrove"),
 ("saltern",         r"saltern|salt ?pan|salt ?field|salina|solar salt|salt crystal|brine pond"),
 ("saline_lake",     r"salt lake|saline lake|hypersaline lake|soda lake|lake"),
 ("marine",          r"marine|sea|ocean|coast|sediment|lagoon|seawater|beach|sand"),
 ("soil",            r"soil|desert|rhizo|plant|root"),
]
def propose(row):
    txt = " ".join(clean(row[c]) for c in ["isolation_source","host","env_broad","env_local","env_medium","organism"]).lower()
    if not txt.strip(): return ("", "no_info")
    for hab, pat in RULES:
        if re.search(pat, txt): return (hab, "keyword")
    return ("", "unmatched")
m[["habitat_proposed","habitat_basis"]] = m.apply(lambda r: pd.Series(propose(r)), axis=1)
m["habitat_final"] = ""      # <- YOU fill/confirm this column
m["coords_source"] = np.where(m.lat.notna(), "biosample", "")

cols = ["accession","strain","pass_qc","organism","country","geo_loc_name","lat","lon","coords_source",
        "isolation_source","host","env_broad","env_local","env_medium","collection_date","biosample",
        "habitat_proposed","habitat_basis","habitat_final"]
m[cols].to_csv(f"{base}/00_metadata/metadata_curation_template.tsv", sep="\t", index=False)

print("Isolates:", len(m), "| pass QC:", int(m.pass_qc.sum()))
print("With coordinates:", int(m.lat.notna().sum()))
print("With country:", int((m.country!="").sum()))
print("Habitat proposed by rule:", int((m.habitat_proposed!="").sum()))
print(m.habitat_proposed.replace("", "UNRESOLVED").value_counts())
print("Countries:", m.country.replace("", "UNRESOLVED").value_counts().head(12).to_dict())
