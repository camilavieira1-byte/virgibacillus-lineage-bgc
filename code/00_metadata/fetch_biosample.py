import os
import urllib.request, urllib.parse, time, re, csv, collections
import xml.etree.ElementTree as ET
import pandas as pd

base = "/PATH/TO/virgibacillus_v2"
m = pd.read_csv(f"{base}/00_metadata/metadata_raw_flagged.tsv", sep="\t")
m = m[m.is_mag == "isolate"]
sam = m[["accession", "biosample"]].dropna()
sam = sam[sam.biosample != "NA"]
acc_by_sam = dict(zip(sam.biosample, sam.accession))
ids = sorted(acc_by_sam)
print("BioSamples to fetch:", len(ids))

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
def fetch(batch):
    q = urllib.parse.urlencode({"db": "biosample", "id": ",".join(batch), "retmode": "xml",
                                "tool": "virgi_revision", "email": os.environ.get("NCBI_EMAIL", "")})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(f"{EUTILS}?{q}", timeout=90) as r:
                return r.read()
        except Exception as e:
            print("retry", attempt, e); time.sleep(3 * (attempt + 1))
    raise RuntimeError("fetch failed")

def dms_to_dec(s):
    s = s.strip()
    # hemisphere + decimal or DMS, e.g. "12.5 N 45.2 E", "12°30'N 45°10'E"
    toks = re.findall(r"(\d+(?:\.\d+)?)\s*(?:°|d)?\s*(?:(\d+(?:\.\d+)?)\s*['′m])?\s*(?:(\d+(?:\.\d+)?)\s*[\"″s])?\s*([NSEWnsew])", s)
    if len(toks) == 2:
        vals = []
        for d, mi, se, h in toks:
            v = float(d) + (float(mi) / 60 if mi else 0) + (float(se) / 3600 if se else 0)
            if h.upper() in "SW": v = -v
            vals.append((h.upper(), v))
        lat = [v for h, v in vals if h in "NS"]; lon = [v for h, v in vals if h in "EW"]
        if lat and lon: return lat[0], lon[0]
    nums = re.findall(r"-?\d+(?:\.\d+)?", s)
    if len(nums) == 2 and not re.search(r"[A-Za-z]", s.replace("NA", "")):
        return float(nums[0]), float(nums[1])
    return None

BAD = {"na", "missing", "not applicable", "not collected", "not provided", "unknown", "none", ""}
rows, attr_names = [], collections.Counter()
for i in range(0, len(ids), 40):
    root = ET.fromstring(fetch(ids[i:i + 40]))
    for bs in root.iter("BioSample"):
        a = bs.get("accession")
        title = (bs.findtext("./Description/Title") or "").strip()
        desc = " ".join((bs.findtext("./Description/Comment/Paragraph") or "").split())
        d = {"biosample": a, "title": title, "comment": desc}
        latlon_raw, lat_sep, lon_sep = "", "", ""
        for at in bs.iter("Attribute"):
            nm = (at.get("attribute_name") or "").lower()
            hn = (at.get("harmonized_name") or "").lower()
            dn = (at.get("display_name") or "").lower()
            val = (at.text or "").strip()
            attr_names[hn or nm] += 1
            if (hn or nm) in ("geo_loc_name", "isolation_source", "collection_date", "host", "env_broad_scale",
                      "env_local_scale", "env_medium", "strain", "sample_type", "type-material", "depth", "temp"):
                d[hn or nm] = val
            if val.lower() in BAD: continue
            if hn == "lat_lon" or re.search(r"lat.*lon", nm + dn):
                latlon_raw = val
            elif "latitude" in nm + dn: lat_sep = val
            elif "longitude" in nm + dn: lon_sep = val
        if not latlon_raw and lat_sep and lon_sep: latlon_raw = f"{lat_sep} {lon_sep}"
        d["latlon_raw"] = latlon_raw
        p = dms_to_dec(latlon_raw) if latlon_raw else None
        d["lat_efetch"], d["lon_efetch"] = (p if p else ("", ""))
        d["latlon_parsed"] = bool(p)
        rows.append(d)
    time.sleep(0.5)

df = pd.DataFrame(rows)
df["accession"] = df.biosample.map(acc_by_sam)
df.to_csv(f"{base}/00_metadata/biosample_efetch.tsv", sep="\t", index=False)
print("Records fetched:", len(df))
print("With a coordinate string:", int((df.latlon_raw != "").sum()))
print("Coordinates parsed to decimal:", int(df.latlon_parsed.sum()))
unparsed = df[(df.latlon_raw != "") & (~df.latlon_parsed)]
print("Coordinate strings NOT parsed (check by hand):", len(unparsed))
print(unparsed[["accession", "latlon_raw"]].to_string(index=False))
print("\nMost common attribute names:", attr_names.most_common(25))
missing = sorted(set(ids) - set(df.biosample))
print("\nBioSamples NOT returned:", len(missing))
for b in missing: print(b, acc_by_sam[b])
if "type-material" in df.columns:
    print("\ntype-material values:"); print(df["type-material"].dropna().value_counts().head(10))
