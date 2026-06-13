# -*- coding: utf-8 -*-
"""GO (BP) + KEGG over-representation analysis on the 43 multi-source drug-disease
intersection targets, via the Enrichr REST API. Reproduces the enrichment described
in Results 3.3 (whose figure/result data were never saved). Writes result CSVs."""
import sys, io, os, json, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import urllib.request, urllib.parse, ssl

ROOT = r"O:\JingTong_CSR_paper"
GENES_CSV = os.path.join(ROOT, "tables", "targets", "intersection_2source.csv")
OUTDIR = os.path.join(ROOT, "tables", "enrichment")
ENRICHR = "https://maayanlab.cloud/Enrichr"
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE

genes = [l.strip() for l in open(GENES_CSV, encoding='utf-8').read().splitlines()[1:] if l.strip()]
print(f"input genes: {len(genes)} -> {','.join(genes[:6])}...")

def add_list(genes):
    boundary = "----EnrichrBoundary7MA4YWxkTrZu0gW"
    body = ""
    for name, val in [("list", "\n".join(genes)), ("description", "JingTong 43 intersection targets")]:
        body += f"--{boundary}\r\nContent-Disposition: form-data; name=\"{name}\"\r\n\r\n{val}\r\n"
    body += f"--{boundary}--\r\n"
    req = urllib.request.Request(ENRICHR + "/addList", data=body.encode("utf-8"),
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    r = json.loads(urllib.request.urlopen(req, timeout=30, context=CTX).read())
    return r["userListId"]

def enrich(uid, library):
    url = f"{ENRICHR}/enrich?userListId={uid}&backgroundType={library}"
    r = json.loads(urllib.request.urlopen(url, timeout=60, context=CTX).read())
    rows = r[library]   # [rank, term, p, zscore, combined, [genes], adjp, oldp, oldadjp]
    out = []
    for x in rows:
        out.append({"rank": x[0], "term": x[1], "p": x[2], "zscore": x[3],
                    "combined_score": x[4], "n_overlap": len(x[5]),
                    "overlap_genes": ";".join(x[5]), "adj_p": x[6]})
    return out

uid = add_list(genes); print("Enrichr userListId:", uid)
time.sleep(1)

def save(rows, name):
    import csv
    path = os.path.join(OUTDIR, name)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote", path, f"({len(rows)} terms)")
    return path

for lib, fn, label in [("KEGG_2021_Human", "KEGG_ORA.csv", "KEGG"),
                       ("GO_Biological_Process_2021", "GO_BP_ORA.csv", "GO-BP"),
                       ("GO_Cellular_Component_2021", "GO_CC_ORA.csv", "GO-CC"),
                       ("GO_Molecular_Function_2021", "GO_MF_ORA.csv", "GO-MF")]:
    rows = enrich(uid, lib); time.sleep(1)
    rows = sorted(rows, key=lambda d: d["adj_p"])
    save(rows, fn)
    print(f"\n=== TOP 12 {label} (by adj p) ===")
    for d in rows[:12]:
        print(f"  adjP={d['adj_p']:.2e}  n={d['n_overlap']:2d}  {d['term']}")
