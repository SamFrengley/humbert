"""
A simple script to check the veracity of Lemma 9.10
"""

from fractions import Fraction as QQ
from sympy import primefactors as prime_divisors
from numpy import prod
from humbert import arithmetic

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
    print("We need to check that the equation on p.43 is >0 for 106 ≤ N ≤ 6000")
    for N in range(106, 6001):
        k = arithmetic.p_adic_val(N, 2)
        M = N // (2 ** k)        
        mu = bounded_mu(M)
        mup = QQ(mu, 2)
        eqn = QQ(N*(N-1)*(N-30), 60) - 3*QQ(2**(2*k), 2)*mup - 6
        assert eqn >= 0
        
    print("Success, now the lemma is proven.")

    
if __name__ == "__main__":
    main()
