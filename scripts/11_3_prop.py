"""
A simple script to check the veracity of Proposition 11.3
"""

from humbert.wnr import WNr_bar

def main():
    """
    main
    """
    rational_cases = [(6,5), (7,3), (8,3), (8,5), (8,7), (9,1), (9,2), (10,1),
                      (10,3), (11,1), (11,3), (12,1), (12,5), (12,7), (12, 11),
                      (13,1), (13,2), (14,1), (14,3), (15,1), (15,2), (15,11),
                      (16,1), (16,3), (16,7), (18,5), (20,11), (24,23)]

    print("We check that: \n  1) (K_Wbar.barC_inf) ≤ -2 for N ≤ 10")
    print("  2) (K_Wbar.barC_inf) ≤ -1 for 11 ≤ N < 24 in the rational cases")
    print("  3) (K_Wbar.barC_inf) = 0 for (N,r)=(24,23)\n")
    for Nr in rational_cases:
        N,r = Nr
        W = WNr_bar(N, r)
        if N <= 10:
            assert W.K_Cinf() <= -2
        elif N >= 11 and N != 24:
            assert W.K_Cinf() <= -1
        elif N == 24:
            assert W.K_Cinf() == 0
    print("Success, now the proposition is proven.")

    
if __name__ == "__main__":
    main()
