import os
import json, time, re, urllib.request, urllib.parse
import pandas as pd
base = "/PATH/TO/virgibacillus_v2/00_metadata"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
EMAIL = os.environ.get("NCBI_EMAIL", "")

def call(tool, **p):
    p.update({"retmode": "json", "tool": "virgi_revision", "email": EMAIL})
    url = E + tool + ".fcgi?" + urllib.parse.urlencode(p)
    for a in range(4):
        try:
            with urllib.request.urlopen(url, timeout=90) as r: 
                time.sleep(0.4); return json.load(r)
        except Exception as e:
            time.sleep(2 * (a + 1))
    return {}

iso = set(pd.read_csv(f"{base}/metadata_raw_flagged.tsv", sep="\t").query("is_mag=='isolate'").accession)
prj = {}
for line in open(f"{base}/summary_all.jsonl"):
    j = json.loads(line); a = j["accession"]
    if a not in iso: continue
    ai = j.get("assembly_info", {})
    p = ai.get("bioproject_accession")
    if not p:
        try: p = ai["bioproject_lineage"][0]["bioprojects"][0]["accession"]
        except Exception: p = "NA"
    prj.setdefault(p, []).append(a)
print("Genomes:", sum(len(v) for v in prj.values()), "| BioProjects:", len(prj))

COORD = re.compile(r"\d+\s*°|\d+\.\d+\s*[NSEW]\b|\blatitude\b|\blongitude\b", re.I)
rows = []
for p, accs in sorted(prj.items(), key=lambda kv: -len(kv[1])):
    d = {"bioproject": p, "n_genomes": len(accs), "accessions": ";".join(accs)}
    if p != "NA":
        uid = (call("esearch", db="bioproject", term=f"{p}[Project Accession]").get("esearchresult", {}).get("idlist") or [None])[0]
        if uid:
            s = call("esummary", db="bioproject", id=uid).get("result", {}).get(uid, {})
            d["title"] = s.get("project_title", ""); d["description"] = (s.get("project_description", "") or "")[:400]
            lk = call("elink", dbfrom="bioproject", db="pubmed", id=uid)
            pm = []
            for ls in lk.get("linksets", []):
                for db in ls.get("linksetdbs", []): pm += db.get("links", [])
            d["pmids"] = ";".join(pm)
            titles, dois = [], []
            if pm:
                ps = call("esummary", db="pubmed", id=",".join(pm[:5])).get("result", {})
                for k in pm[:5]:
                    x = ps.get(k, {}); titles.append(x.get("title", ""))
                    dois += [a["value"] for a in x.get("articleids", []) if a.get("idtype") == "doi"]
            d["paper_titles"] = " | ".join(titles); d["dois"] = ";".join(dois)
            d["coords_in_text"] = bool(COORD.search(d["title"] + " " + d["description"]))
    rows.append(d)

df = pd.DataFrame(rows)
df.to_csv(f"{base}/bioproject_papers.tsv", sep="\t", index=False)
print("Projects with a linked PubMed paper:", int((df.get("pmids", "") != "").sum()))
print("Genomes covered by a project with a paper:", int(df.loc[df.get("pmids", "") != "", "n_genomes"].sum()))
print("Projects with coordinate-like text:", int(df.get("coords_in_text", False).sum()))
print(df[["bioproject", "n_genomes", "title", "dois"]].head(25).to_string(index=False, max_colwidth=60))
