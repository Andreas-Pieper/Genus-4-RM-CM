from lmfdb import db

def read_params(t):
    D = []
    if t==7:
        file_in = '/home/sschiavo/github/Genus-4-RM-CM/examples/Shimura_237_param.m'
    elif t==9:
        file_in = '/home/sschiavo/github/Genus-4-RM-CM/examples/Shimura_239_param.m'
    else:
        raise ValueError('t must be 7 or 9')
    with open(file_in,'r') as f_in:
        lines = f_in.readlines()
        for line in lines:
            spl = line.split('|')
            if len(spl) > 3: # a couple entries have no parameters
                print("Reading %s" % spl[1])
                D.append((spl[1],spl[2]))
    return list(set(D))

def sort_fields(t,D):
    D_new = []
    for pair in D:
        rec = db.nf_fields.lookup(pair[0])
        assert rec['disc_abs'] % t == 0
        d = (rec['disc_sign']*rec['disc_abs']).factor()
        D_new.append((pair[0],d,pair[1]))
        # sort by discriminant
        D_new = sorted(D_new, key = lambda el: el[1])
    return D_new

# pair consists of an LMFDB field label and a stringified polynomial; t is either 7 or 9, corresponding to the (2,3,7) or (2,3,9) triangle group
def table_row(trip,t):
    fld_lab, d, inv_poly = trip 
    print("Writing row for %s" % fld_lab)
    # polredabs polynomial
    R.<x> = PolynomialRing(QQ)
    f = sage_eval(inv_poly, locals = locals())
    fabs = R(str(gp.polredabs(gp.polredbest(f))))
    cs = [ZZ(el) for el in fabs.list()]
    inv_lab = db.nf_fields.lucky({'coeffs':cs}, projection='label')
    # maybe replace poly with its LMFDB field label
    # currently there are two deg 15 polys whose numflds aren't in LMFDB
    link = "https://www.lmfdb.org/NumberField/%s" % fld_lab
    if inv_lab:
        link_inv = "https://www.lmfdb.org/NumberField/%s" % inv_lab
        s = r"\href{%s}{%s} & $%s$ & \href{%s}{%s}\\" % (link, fld_lab, latex(d), link_inv, inv_lab)
    else:
        s = r"\href{%s}{%s} & $%s$ & $%s$\\" % (link, fld_lab, latex(d), latex(fabs))
    return s

def make_table(t, D, file_out):
    R.<x> = PolynomialRing(QQ)
    with open(file_out,'a') as f_out:
        # table header
        f_out.write(r'LMFDB label & $\disc(K)$ & Field of definition\\'+'\n')
        for trip in D:
            print("Processing %s" % trip[0])
            f_out.write(table_row(trip,t)+"\n")
    return "Table written to %s" % file_out

def missing_fields(t, file_out):
    flds = read_params(t)
    flds = sort_fields(t,flds)
    with open(file_out,'a') as f_out:
        for trip in flds:
            fld_lab, d, inv_poly = trip 
            print("Checking %s" % fld_lab)
            # polredabs polynomial
            R.<x> = PolynomialRing(QQ)
            f = sage_eval(inv_poly, locals = locals())
            fabs = R(str(gp.polredabs(gp.polredbest(f))))
            cs = [ZZ(el) for el in fabs.list()]
            inv_lab = db.nf_fields.lucky({'coeffs':cs}, projection='label')
            if not inv_lab:
                print("%s missing" % fabs)
                f_out.write(str(fabs)+"\n")
    return "Missing polys written to %s" % file_out


# now do it
D7 = read_params(7)
D7 = sort_fields(7,D7)
D9 = read_params(9)
D9 = sort_fields(9,D9)
print("Making table for 2,3,7")
make_table(7,D7,"table7.txt")
print("Making table for 2,3,9")
make_table(9,D9,"table9.txt")
#missing_fields(7,"missing7.txt")
#missing_fields(9,"missing9.txt")
