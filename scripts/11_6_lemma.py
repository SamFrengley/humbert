"""
A simple script to check the veracity of Lemma 11.2
"""

from humbert import arithmetic
from humbert.wnr import WNr, WNr_o

def main():
    """
    Main
    """
    print("\nWe wish to show that K_small^2 > 0 for those (N,r) where")
    print("  p_g(W_(N,r)) ≥ 3 and N < 106 excepting the cases where")
    print("  (N,r) = (24,1), (24,5), (24,7), (24,17). In these cases we must")
    print("  show that K_small^2 = 0.\n")

    exceptions = [(24,1), (24,5), (24,7), (24,17)]
    
    for N in range(20,106):
        for r in arithmetic.possible_r(N):
            W = WNr(N, r)
            if W.p_g() >= 3:
                Wo = WNr_o(N,r)
                flag1 = ((N, r) in exceptions)
                flag2 = (Wo.K_K() > 0)
                assert flag1 or flag2

    for Nr in exceptions:
        W = WNr_o(Nr[0], Nr[1])
        assert W.K_K() == 0
                
    print("Success, now the lemma is proven.\n")

    
if __name__ == "__main__":
    main()
