from lmfdb import db

# fld_lab is an LMFDB field label; t is either 7 or 9, corresponding to the (2,3,7) or (2,3,9) triangle group
def table_row(fld_lab,t):
    rec = db.nf_fields.lookup(fld_lab)
    if rec['disc_abs'] % t == 0:
        d = (rec['disc_sign']*rec['disc_abs']).factor()
        link = "https://www.lmfdb.org/NumberField/%s" % fld_lab
        if rec['galois_label'] == '8T13':
            gal_str = r'C_2 \times A_4'
        if gal_str:
            s = r"\href{%s}{%s} & $%s$ & $%s$\\" % (link, fld_lab, latex(d), gal_str)
            return s
    return ''

def make_table(t, file_out, file_in='/home/sschiavo/github/Genus-4-RM-CM/CM/curves_class_8t13.txt'):
    R.<x> = PolynomialRing(QQ)
    #R = PolynomialRing(QQ, 'x')
    #x = R.gens()[0]
    with open(file_out,'a') as f_out:
        # table header
        f_out.write(r'LMFDB label & $\disc(K)$ & $\Gal(K/\Q)$\\'+'\n')
        # look for entries corresponding to CM points on Mumford Shimura curve
        with open(file_in,'r') as f_in:
            lines = f_in.readlines()
            lines = list(set(lines)) # eliminate redundancies
            for line in lines:
                line = line.strip()
                spl = line.split('|')
                fld_lab = spl[0]
                print("Processing %s" % fld_lab)
                crv_type = spl[1]
                # Check for QQ(zeta_7)^+ or QQ(zeta_9)^+ being contained in reflex field 
                rec = db.nf_fields.lookup(fld_lab)
                # reflex field should have degree 6
                rfs = [el for el in db.nf_fields_reflex.search({'nf_label':fld_lab}) if len(el['rf_coeffs'])==7]
                if (t == 7) and (crv_type == 'hyp'):
                    g = x^3 - x^2 - 2*x + 1 # defining poly of QQ(zeta_7)^+
                elif (t == 9) and (crv_type == 'true'):
                    g = x^3 - 3*x - 1 # defining poly of QQ(zeta_9)^+
                else:
                    continue
                zeta_bool = False
                for rf in rfs:
                    F.<a> = NumberField(R(rf['rf_coeffs']))
                    #F = NumberField(R(rf['rf_coeffs']),'a')
                    #a = F.gens()[0]
                    if len(g.roots(F)) != 0:
                        zeta_bool = True
                if zeta_bool:
                    f_out.write(table_row(fld_lab,t)+"\n")
    return "Table written to %s" % file_out
