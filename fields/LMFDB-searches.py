# QQ(zeta_3) has label 2.0.3.1
# QQ(zeta_3) has coeffs 1.-1.1
# QQ(zeta_5) has label 4.0.125.1
# QQ(zeta_5) has coeffs 1.-1.1.-1.1

non_dJN = db.nf_fields.search({'degree':8, 'cm':True, 'subfields':{'$notcontains':['1.-1.1', '1.-1.1.-1.1']}})
dJN = db.nf_fields.search({'degree':8, 'cm':True, 'subfields':{'$contains':['1.-1.1']}})
dJN = list(dJN)
dJN2 = db.nf_fields.search({'degree':8, 'cm':True, 'subfields':{'$contains':['1.-1.1.-1.1']}})
dJN2 = list(dJN2)
dJN = dJN + dJN2

with open("/scratch/home/sschiavo/github/Genus-4-RM-CM/CM/deJongNoot-fields.txt", "a") as my_file:
    for rec in dJN:
        s = "%s|%s" % (rec['label'], rec['coeffs'])
        my_file.write(s.replace(" ",""))
