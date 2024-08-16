"""
This is a script which outputs the tables in Appendix B.
"""

from humbert.znr import *
from humbert.wnr import *

def table_str(N):
    """
    Get a string which tex outputs as one line of the table
    """
    poss_r = possible_r(N)
    all_r = []
    for r in poss_r:
        if Ztil(N, r).kodaira_dim() != -1:
            all_r.append(r)
            
    num_r = len(all_r)
    ret = f"    \\multirow{{{num_r}}}{{*}}{{${N}$}} "
    for r in all_r:
        W = WNr(N, r)
        Wbar = WNr_bar(N, r)
        Z = Ztil(N, r)
        pg_Z = Z.p_g()
        kd_Z = Z.kodaira_dim()
        pg_W = W.p_g()
        kd_W = W.kodaira_dim()
        K_Cinf = Wbar.K_Cinf()
        KK_W = W.K_K()
        
        if kd_W != -1:
            Wo = WNr_o(N, r)
            KK = Wo.K_K()
        else:
            KK = " "

        if rho(N, -r) == 0:
            line = f"& ${r}$ & ${pg_Z}$ & ${kd_Z}$ & ${pg_W}$ "
            line += f"& ${K_Cinf}$ & ${KK_W}$ & ${KK}$ & ${kd_W}$"
            line += "\\\\ \n"
        else:
            line = f"& $\\boldsymbol{{{r}}}$ & ${pg_Z}$ & ${kd_Z}$ & ${pg_W}$ "
            line += f"& ${K_Cinf}$ & ${KK_W}$ & ${KK}$ & ${kd_W}$"
            line += "\\\\ \n"
            
        if r == 1:
            ret += line
        else:
            ret += "    "
            ret += line
            
    return ret
    
def generate_table(low, high):
    """
    Get the string for the entire table
    """
    caption = "Numerical invariants of the surfaces $\\WNr{N}{r}$. The \n"
    caption += "    surface $\\WNro{N}{r}$ is defined only when $\\WNr{N}{r}$ is \n"
    caption += "    not rational, so the column recording $\\KWo^2$ is left \n"
    caption += "    empty if $\\WNr{N}{r}$ is rational. For the bolded values of \n"
    caption += "    $r$ the surface $\\WNr{N}{r}$ is birational to the Humbert\n"
    caption += "    surface $\\mathcal{H}_{N^2}$ in Theorem 1.1."
    
    ret = "\\begin{table}[H] \n"
    ret += "  \\centering \n"
    ret += "  \\begin{tabular}{cc|ccccccc} \n"
    ret += "    $N$ & $r$ & $p_g(\\ZNrtil{N}{r})$ & $\\kappa(\\ZNrtil{N}{r})$ "
    ret += "& $p_g(\\WNr{N}{r})$ & $K_{\\overbar{W}} \\cdot \\overbar{C}_\\infty$ "
    ret += "& $K_W^2$ & $\\KWo^2$ & $\\kappa(\\WNr{N}{r})$ \\\\ \n"
    ret += "    \\hline\n"

    for N in range(low, high):
        ret += table_str(N)
        if N != high - 1:
            ret += "    \\hdashline\n"
        
    ret += "  \\end{tabular} \n"
    ret += f"  \\caption{{{caption}}} \n"
    ret += "  \\label{table:placeholder} \n"
    ret += "\\end{table}"
    return ret

def main():
    """
    main
    """
    print(generate_table(6, 21))
    print("\n\n")
    print(generate_table(21, 34))
    
if __name__ == "__main__":
    main()
