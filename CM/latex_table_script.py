from lmfdb import db

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
    with open(file_out,'a') as f_out:
        # table header
        f_out.write(r'LMFDB label & $\disc(K)$ & $\Gal(K/\Q)$\\'+'\n')
        with open(file_in,'r') as f_in:
            lines = f_in.readlines()
            lines = list(set(lines))
            for line in lines:
                line = line.strip()
                spl = line.split('|')
                fld_lab = spl[0]
                print("Processing %s" % fld_lab)
                crv_type = spl[1]
                if t == 7:
                    if crv_type == 'hyp':
                        f_out.write(table_row(fld_lab,t)+"\n")
                if t == 9:
                    if crv_type == 'true':
                        f_out.write(table_row(fld_lab,t)+"\n")
    return "Table written to %s" % file_out
