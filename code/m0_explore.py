
import sys, json, time
sys.path.insert(0, __import__("os").path.dirname(__file__))
import lib
api = lib.AikiSolAPI(delay=1.0)
ex = "MTKYEESMKAKSGVGPKEYIKQLNDSFVKAQSASGVEVVTGATHSSESFQNYAQQLIQAAQAGNTDTIEIDNGATLKDGTYSLKEKNDSNGYHTTFSMTVKDGKVTESNYDNVNADGKSKKDDTEYESKMKDVAGVGPKEYIETLNKEFVKAMGEEDGSPAGVEVVTGATHSTHSFINYAQQLVNAAEKGDTTEIVVDNIVTK"
for name, seq in [("landing-page example", ex), ("doc example 20-mer", "MKLITVLVLALLAVAVAFPV")]:
    st, tx, err = api.lookup_raw(seq, tag="M0_format")
    print("== lookup", name, st, err); print(tx[:1500])
fa = ">p1\n"+ex+"\n>p2\nMKLITVLVLALLAVAVAFPV\n"
for fmt in ("json","csv"):
    st, tx, err = api.predict_fasta_raw(fa, fmt, tag="M0_format")
    print("== predict_fasta", fmt, st, err); print(tx[:1500])
print("requests so far:", lib.api_count())
