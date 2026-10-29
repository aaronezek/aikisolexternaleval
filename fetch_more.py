import requests, json, os
D = "data/"
r = requests.get("https://raw.githubusercontent.com/Violet969/PLM_Sol/main/embedding_dataset/Train_dataset.fasta", timeout=120); r.raise_for_status()
open(D+"plmsol_train_dataset.fasta", "wb").write(r.content); print("train fasta", len(r.content))
j = json.load(open(D+"zenodo_record.json")); base = {f["key"]: f["links"]["self"] for f in j["files"]}
k = "apache_tier__aikisol_canonical_147k_train.csv"
r = requests.get(base[k], timeout=600); r.raise_for_status(); open(D+"zenodo/"+k, "wb").write(r.content); print("aikisol train csv", len(r.content))
