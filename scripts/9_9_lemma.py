"""
A simple script to check the veracity of Lemma 9.9
"""

from fractions import Fraction as QQ
from sympy import primefactors as prime_divisors
from numpy import prod
from humbert import arithmetic
from humbert import wnr

def bounded_mu(M):
    """
    Since mu = M**2 * prod(( p + (-r/p) ) / p) we have
    mu <= M**2 * prod((p + 1)/p) and this is what we return
    """
    if M == 1:
        return 1
    
    mu = M**2 * prod(
        [
            QQ(int(p + 1), int(p))
            for p in prime_divisors(M)
        ]
    )
    return mu


def main():
    """
    Main
    """
    print("We first need that the equation on p.43 is ≥3 for 67 ≤ N ≤ 2500")
    for N in range(67, 2501):
        k = arithmetic.p_adic_val(N, 2)
        M = N // (2 ** k)
        mu = bounded_mu(M)
        mup = QQ(mu, 2)
        eqn = QQ(N*(N-1)*(N-23), 480) - QQ(2**(2*k), 8)*mup - 1
        eqn = QQ(eqn, 2)
        assert eqn >= 3
    print("Success, the above didn't throw an error.")
    
    print("We need to check p_g(W_(N,r)) ≥3 for N ≤ 67 not in Thm 1.10 (i)--(iii).")
    # note that the exceptions in Thm 1.10 are all those N < 20, together
    # with the following
    exceptions = [(20, 1), (20, 3), (20, 11), (21, 2), (21, 5), (22, 1),
                  (24, 11), (24, 23)]

    for N in range(20, 67):
        for r in arithmetic.possible_r(N):
            W = wnr.WNr(N, r)
            flag1 = ((N, r) in exceptions)
            flag2 = (W.p_g() >= 3)
            assert flag1 or flag2

    print("Success, now the lemma is proven.")

    
if __name__ == "__main__":
    main()
