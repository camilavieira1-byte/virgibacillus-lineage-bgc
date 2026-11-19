import urllib.request, json, time, re
import pandas as pd
base = "/PATH/TO/virgibacillus_v2/00_metadata"
m = pd.read_csv(f"{base}/metadata_raw_flagged.tsv", sep="\t")
m = m[m.is_mag == "isolate"][["accession", "biosample"]]
got = set(pd.read_csv(f"{base}/biosample_efetch.tsv", sep="\t").biosample)
todo = m[~m.biosample.isin(got)]
print("BioSamples to query at EBI:", len(todo))

BAD = {"na", "missing", "not applicable", "not collected", "not provided", "unknown", "none", ""}
rows = []
for r in todo.itertuples():
    url = f"https://www.ebi.ac.uk/biosamples/samples/{r.biosample}"
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=60) as f:
            j = json.load(f)
    except Exception as e:
        rows.append({"accession": r.accession, "biosample": r.biosample, "error": str(e)}); time.sleep(0.5); continue
    ch = {k.lower(): (v[0].get("text", "") if v else "") for k, v in j.get("characteristics", {}).items()}
    d = {"accession": r.accession, "biosample": r.biosample, "ebi_name": j.get("name", ""), "error": ""}
    for k, v in ch.items():
        if re.search(r"\blat|\blon|coordinate", k): d["coord_" + k] = v
    for k in ("geographic location (country and/or sea)", "geo_loc_name", "isolation source host-associated",
              "isolation_source", "isolation source", "collection date", "strain", "host", "sample type", "description", "title"):
        if k in ch: d[k] = ch[k]
    rows.append(d); time.sleep(0.4)

df = pd.DataFrame(rows)
df.to_csv(f"{base}/biosample_ebi.tsv", sep="\t", index=False)
coordcols = [c for c in df.columns if c.startswith("coord_")]
has = df[coordcols].apply(lambda s: s.fillna("").astype(str).str.lower().isin(BAD) == False).any(axis=1) if coordcols else pd.Series(False, index=df.index)
print("Queried:", len(df), "| errors:", int((df.get("error", "") != "").sum()))
print("With some coordinate field:", int(has.sum()))
print("Coordinate-like fields found:", coordcols)
print(df[["accession", "biosample"] + coordcols].head(30).to_string(index=False))
