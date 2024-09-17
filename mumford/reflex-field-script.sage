from lmfdb import db

load("/scratch/home/sschiavo/github/Genus-4-RM-CM/mumford/mumford_jac.py")
R.<x> = PolynomialRing(QQ)
for lab in CM_flds:
  print(f"Computing for label {lab}")
  refs = db.nf_fields_reflex.search({'nf_label': lab})
  for ref in refs:
    cs = ref['rf_coeffs']
    print(f"cs length = {len(cs)}")
    if len(cs) != 7:
      continue
    K.<nu> = NumberField(R(cs))
    G = K.galois_group()
    print(G.transitive_label())
    if G.transitive_label() == '6T11':
      CM_cs = db.nf_fields.lookup(lab, projection='coeffs')
      s = "|".join([str(el) for el in [lab, CM_cs, cs]])
      s = s.replace(" ","")+"\n"
      print(s)
      with open("/scratch/home/sschiavo/github/Genus-4-RM-CM/mumford/mumford-results.txt",'a') as f:
        f.write(s)