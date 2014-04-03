def get_pdflabel(s):
    pdf_labels = {}
    if s in pdf_labels:
        return pdf_labels[s]
    elif 'HERAMCDR_DISCMSJETSCMSWASYM' in s:
        return 'HERA DIS + CMS WASYM/Jets'
    elif 'HERAMCDR_DISCMSJETS' in s:
        return 'HERA DIS + CMS Jets'
    elif 'HERAMCDR_DISCMSWASYM' in s:
        return 'HERA DIS + CMS W ASYM'
    elif 'HERAMCDR_DIS' in s:
        return 'HERA DIS'
    elif 'HERADISCMSJETS' in s:
        return 'HERA DIS + CMS Jets'
    elif 'HERADISATLASJETS' in s:
        return 'HERA DIS + Atlas Jets'
    elif 'HERADISCMSATLASJETS' in s:
        return 'HERA DIS + Atlas + CMS Jets'
    elif 'HERADIS' in s:
        return 'HERA DIS'
    else:
        return s.replace('_','\_')


def get_partonlabel(s, short=False):
    parton_labels = {0: "gluon",
                     1: "d quark",
                     2: "u quark",
                     3: "s quark",
                     4: "c quark",
                     -1: 'd antiquark',
                     -2: 'u antiquark',
                     7: '$d_{\mathrm{val}}$ quark',
                     8: '$u_{\mathrm{val}}$ quark',
                     9: 'sea quarks'}

    parton_labels_short = {0: "$g$",
                           1: r"$D$",
                           2: r"$U$",
                           3: r"$s$",
                           4: r"$c$",
                           -1: r'$\bar D$',
                           -2: r'$\bar U$',
                           7: r'$d_{V}$',
                           8: r'$u_{V}$',
                           9: r'$\mathrm{sea}$'}

    if short:
        if s in parton_labels:
            return parton_labels_short[s]
        else:
            return s
    else:
        if s in parton_labels:
            return parton_labels[s]
        else:
            return s


def get_q2label(q2):
    q2_str = r'{:.2G}'.format(q2)
    if "E" in q2_str:
        base, exponent = q2_str.split('E')
        return r"{0}\times 10^{{{1}}}".format(base, int(exponent))
    else:
        return q2_str

def get_plot_scalefactor(flavor, q2=1.9):

    scalefactors = {
            0 : { 1.9 : 0.2},
            1  : { 1.9 : 1.0},
            2  : { 1.9 : 1.0},
            7  : { 1.9 : 1.0},
            8  : { 1.9 : 1.0},
            9  : { 1.9 : 0.2},
            }

    return scalefactors[flavor].get(q2, 1.0)
