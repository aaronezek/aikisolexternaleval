
import requests, json, hashlib, time
j=json.load(open("aikisol_eval/data/zenodo_record.json"))
f=[x for x in j["files"] if x["key"]=="apache_tier__aikisol_v2_canonical_147k_full.pt"][0]
url=f["links"]["self"]; out="aikisol_eval/data/zenodo/aikisol_v2_canonical_147k_full.pt"
print(f["checksum"], f["size"], flush=True)
t=time.time(); h=hashlib.md5(); n=0
with requests.get(url, stream=True, timeout=120) as r:
    r.raise_for_status()
    with open(out,"wb") as fh:
        for ch in r.iter_content(1<<22):
            fh.write(ch); h.update(ch); n+=len(ch)
            if n % (1<<28) < (1<<22): print(n/1e6,"MB", round(time.time()-t), "s", flush=True)
print("done", n, "md5", h.hexdigest(), "expected", f["checksum"], flush=True)
