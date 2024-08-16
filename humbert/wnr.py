"""
wnr.py
======

This module provides functionality to compute various invariants associated to
surfaces in the birational equivalence class of the surface W_(N,r) -- which is
itself birational to Z_(N,r)^sym

Classes
-------
- WNrCusp(N, r, d, q)
- WNr(N, r)
- WNr_bar(N, r)
- WNr_o(N, r)

Examples
--------
Some basic usage:

>>> from humbert.wnr import WNr, WNr_bar, WNr_o
>>> W = WNr(12, 1)
>>> W.p_g()
0
>>> W.K_Fg()                                                # [ K.(F_(c,w)^+)* ]
[0]

>>> W = WNr_bar(24, 23)
>>> W.p_g()
0
>>> W.K_Cinf()                                              # nearly win already
0
>>> W.kodaira_dim()
-1

>>> W = WNr_o(22, 7)
>>> W.kodaira_dim()
2
>>> W.K_K()                                                 # > 1 so we win
1

>>> W = WNr_o(19, 1)                                        # W^small_(19,1)
>>> W.kodaira_dim()                                         # prop. ell
1
>>> print(W)
The surface W^small_(19,1)
κ     = 1
p_g   = 2
K^2   = -1

"""

################################################################################
# External modules
from math import gcd
from sympy import totient as euler_phi, divisors

################################################################################
# Local modules
from .arithmetic import *
from .classnumbers import class_number                      # hardcoded h(-d)
from .modularcurves import *
from .znr import *

################################################################################
# This module starts here

class WNrCusp:
    """Cusp E_(∞,d,q)^* on W_(N,r).

    Parameters
    ----------
    N : int
        Level of HMS.
    r : int
        Power of the congruence.
    d : int
        Integer dividing `N`, singularity of type (`d`,`q`).
    q : int
        Integer coprime to `d` and so that `q`*`r` is a square mod d.

    Arguments
    ---------
    N : int
        Level of HMS.
    r : int
        Power of the congruence.
    d : int
        Integer dividing `N`, singularity of type (`d`,`q`).
    q : int
        Integer coprime to `d` and so that `q`*`r` is a square mod d.
    Z_cusp : ZtilCusp
        The corresponding cusp on ~Z_(N,r) (more precisely, one of its connected
        components in the analytic topology).

    """

    def __init__(self, N, r, d, q):
        self.N = N
        self.r = r
        self.d = d
        self.q = q
        self.Z_cusp = ZtilCusp(N, r, d, q)
        
    def self_ints(self):
        """Compute the components and self-intersections E_(∞,d,q)^* on W_(N,r).

        Returns
        -------
        list 
            The self intersections of the components of E_(∞,d,q)^*
        
        """
        if not self.Z_cusp.is_fixed:
            return self.Z_cusp.ctd_frac
        else:
            err_str = "Sorry, I have not implemented this. Email me if you\n"
            err_str += "need it. The reason is that the resolution meets a\n"
            err_str += "fixed curve, so the self-intersections change."
            raise NotImplementedError(err_str)

class WNr:
    """The class of the surface W_(N,r)

    Parameters
    ----------
    N : int
        Level of the Humbert surface
    r : int
        Power of the congruence

    Attributes
    ----------
    N : int
        Level of the Humbert surface
    r : int
        Power of the congruence
    M : int
        Odd part of N
    k : int
        So that N = M * 2**k
    
    """

    def __init__(self, N, r):
        if gcd(N, r) != 1:
            raise ValueError("N,r should be coprime")
        if Ztil(N, r).kodaira_dim() == -1:
            raise ValueError(
                "W_(N,r) is rational but not defined since Z_(N,r) is rational"
            )
        self.N = N
        self.r = r
        self.k = p_adic_val(N, 2)
        self.M = N // 2 ** (self.k)

    def __str__(self):
        ret = "The surface W_({N},{r})\n".format(N=self.N, r=self.r)
        ret += "κ     = {}\n".format(self.kodaira_dim())
        ret += "p_g   = {}\n".format(self.p_g())
        ret += "K^2   = {}\n".format(self.K_K())
        # ret += "χ_top = {}".format(self.chi_top())
        return ret

    def p_g(self):
        """BIRATIONAL INVARIANT: The geometric genus of W_(N,r).

        Returns
        -------
        int
            The geometric genus of W_(N,r).

        Notes
        -----
        This is just an implementation of the formula in [1]_ Theorem 9.3 (cf.
        [2]_ Equation (10)).

        References
        ----------
        .. [1] S. Frengley, 'On the geometry of the Humbert surface of
           discriminant N^2', this article (2024).
        .. [2] C.F. Hermann, 'Symmetrische Modulflachen der Diskriminante p^2',
           Manuscripta Math. 76 (2992), no.2, 169-179.
        
        """
        N = self.N
        r = self.r
        Z = Ztil(N, r)
        pg_Z = Z.p_g()
        k = self.k
        M = self.M

        X = X_w(M, r)
        mu = X.mu
        e_inf = X.e_inf

        if pg_Z == 0:
            ret = 0
        else:
            if k == 0:
                ret = 24 * pg_Z
                ret += 9 * rho(N, r)
                ret += 6 * rho(N, 2 * r)
                ret += 4 * rho(N, 3 * r)
                ret -= mu
                ret += 6 * e_inf
                ret -= 24
                ret = ret // 48
            elif k == 1:
                ret = 12 *pg_Z
                ret += 6 * rho(M, r) 
                ret += 2 * rho(M, 3 * r) 
                ret -= 2 * mu 
                ret += 6 * e_inf
                ret -= 12 
                ret = ret // 24
            elif k == 2:
                if r % 4 == 1:
                    ret = 4 * pg_Z
                    ret += 3 * rho(M, r) 
                    ret += 3 * e_inf 
                    ret -= 2 * mu
                    ret -= 4
                    ret = ret // 8
                elif r % 4 == 3:
                    ret = 12 * pg_Z
                    ret += + 3 * rho(M, r)
                    ret += 4 * rho(M, 3 * r)
                    ret += 15 * (e_inf)
                    ret -= 10 * mu
                    ret -= 12
                    ret = ret // 24
            elif k >= 3:
                if r % 8 == 1:
                    if k == 3:
                        const = 6
                    else:
                        const = 8
                    ret = 8 * pg_Z
                    ret += 2 * const * rho(M, r)
                    ret += 3 * 2 ** (k - 1) * e_inf
                    ret -= 2 ** (2 * k - 2) * mu
                    ret -= 8
                    ret = ret // 16
                elif r % 8 == 3:
                    ret = 24 * pg_Z
                    ret += 2 * 8 * rho(M, 3 * r)
                    ret += 3 * 2 ** (k + 1) * e_inf
                    ret -= 2 ** (2 * k) * mu
                    ret -= 24
                    ret = ret // 48
                elif r % 8 == 5:
                    if k == 3:
                        const = 2
                    else:
                        const = 0
                    ret = 8 * pg_Z
                    ret += 2 * const * rho(M, r)
                    ret += 3 * 2 ** (k - 1) * e_inf
                    ret -= 2 ** (2 * k - 2) * mu
                    ret -= 8
                    ret = ret // 16
                elif r % 8 == 7:
                    ret = 8 * pg_Z 
                    ret += 3 * 2 ** (k) * e_inf 
                    ret -= 2 ** (2 * k - 1) * mu
                    ret -= 8
                    ret = ret // 16
        return ret

    def K_K(self):
        """Self intersection of the canonical on W_(N,r)

        Returns
        -------
        int
            The self intersection, K.K, of a canonical divisor on W_(N,r)

        Notes
        -----
        This is just an implementation of the formula in [1]_ Lemma 9.4.

        References
        ----------
        .. [1] S. Frengley, 'On the geometry of the Humbert surface of
           discriminant N^2', this article (2024).       
        
        """
        N = self.N
        r = self.r
        M = self.M
        k = self.k
        if k == 0:
            ret = 13 * rho(M, r) 
            ret += 2 * 4 * rho(M, 2 * r) 
            ret += 2 * 2 * rho(M, 3 * r)
            ret = ret // 2
        elif k == 1:
            ret = 9 * rho(M, r) + 2 * rho(M, 3 * r)
        elif k == 2 and r % 4 == 1:
            ret = 14 * rho(M, r)
        elif k == 2 and r % 4 == 3:
            ret = 4 * rho(M, r) + 4 * rho(M, 3 * r)
        elif k == 3 and r % 8 == 1:
            ret = 28 * rho(M, r)
        elif k >= 3 and r % 8 == 3:
            ret = 8 * rho(M, 3 * r)
        elif k == 3 and r % 8 == 5:
            ret = 8 * rho(M, r)
        elif k >= 4 and r % 8 == 1:
            ret = 36 * rho(M, r)
        else:
            ret = 0
        # at this point 2Kw^2 - Kz^2 + 2KzF - F^2 = ret
        # => Kw^2 = 1/2*Kz^2 - 3/2*KzF + sum(g - 1) + ret/2
        Z = Ztil(N, r)
        X = Xg_plus(N, r)
        genera = X.genus()
        gg = sum([g - 1 for g in genera])
        KZ2 = Z.K_K()
        KzF = sum(Z.K_Fg())
        return (KZ2 - 3 * KzF + 2 * gg + ret) // 2

    def kodaira_dim(self):
        """BIRATIONAL INVARIANT: The Kodaira dimension of W_(N,r).

        Returns
        -------
        int
            The Kodaira dimension of W_(N,r).
        
        """
        return min([2, self.p_g() - 1])

    def chi_top(self):
        """The topological Euler characteristic of W_(N,r).

        Returns
        -------
        int
            The integer χ_top(W_(N,r)).

        Raises
        ------
        NotImplementedError
            Because I haven't implemented this yet. Please contact me if you
            need it.
        
        """
        err_str = "Sorry not implemented. This should follow from Max \n"
        err_str += "Noether's formula. If you need this please email me."
        raise NotImplementedError(err_str)           
    
    def K_Cinf(self):
        """The intersection K.C_∞* of C_∞* with a canonical divisor on WNr.

        Returns
        -------
        int
            The intersection number K.C_∞*.

        Notes
        -----
        This is just an implementation of the formula in [1]_  Lemma 11.1.

        References
        ----------
        .. [1] S. Frengley, 'On the geometry of the Humbert surface of
           discriminant N^2', this article (2024).       
        
        """
        N = self.N
        r = self.r
        Z = Ztil(N, r)
        ret = Z.K_Cinf()
        ret += -euler_phi(N) // 2
        return ret

    def K_Fm(self, m):
        """An upper bound on K.F_(m,λ)*.

        Parmaeters
        ----------
        m : int
            An integer, the level of the HZ divisor F_(m,λ)^*

        Returns
        -------
        ret : int
            (A bound on) the intersection number K.F_m*
        flag : bool
            A boolean which is True if F_m* is a nonsingular curve and
            K.F_m* = ret

        Notes
        -----
        This is [1]_ Lemma 10.1.

        References
        ----------
        .. [1] S. Frengley, 'On the geometry of the Humbert surface of
           discriminant N^2', this article (2024).        
        
        """
        N = self.N
        r = self.r
        X = X_0(m)
        mu = X.mu
        nu_3 = X.e_3
        nu_inf = X.e_inf
        if m % 4 == 3:
            Fm_FF = class_number(-4 * m) + class_number(-m)
        else:
            Fm_FF = class_number(-4 * m)
        ret = (Ztil(N, r).K_Fm(m) - Fm_FF) // 2
        if m < N:
            flag = True
        elif ret <= -1:
            if self.p_g() > 0:
                flag = True
            else:
                flag = False
        else:
            flag = False
        return (ret, flag)

    def K_Fg(self):
        """The intersection number K.F_g* of F_g* with a canonical divisr.

        More precisely this gives the K.F_g where F_g ranges over those F_g
        which were fixed on Z_(N,r).

        Returns
        -------
        list
            The intersection numbers [K.F_g*] for each component of FF. The
            order of the components is determined by the following.
            - `self.k` == 0 : [w]
            - `self.k` == 1 : [(I,w), (c,w)]
            - `self.k` == 2 and r == 1 : [(c,w)]
                            and r == 3 : [(ns, w), (s, w), (c,w)]
            - `self.k` >= 3 and r == 1 : [(c,w)]
            - `self.k` >= 3 and r == 3 : [(ns, w), (ωns, w), (c,w)]
            - `self.k` >= 3 and r == 5 : [(c,w)]
            - `self.k` >= 3 and r == 7 : [(s, w), (ωs, w), (c,w)]

        Notes
        -----
        This is just from the projection formula.
        
        """
        N = self.N
        r = self.r
        k = self.k
        M = self.M
        Z = Ztil(N, r)
        X = Xg_plus(N, r)
        Ktil = Z.K_Fg()
        genera = X.genus()
        # the following are giving K_Ztil.~F_g - KZo.Fg_o
        if k == 0:
            ret = [(3 * rho(M, r) + 2 * rho(M, 2 * r) + rho(M, 3 * r)) // 2]
        elif k == 1:
            ret = [(rho(M, 3 * r) + 3 * rho(M, r)) // 2, rho(M, r) // 2]
        elif k == 2 and r % 4 == 1:
            ret = [3 * rho(M, r)]
        elif k == 2 and r % 4 == 3:
            ret = [rho(M, 3 * r), 0, rho(M, r)]
        elif k == 3 and r % 8 == 1:
            ret = [6 * rho(M, r)]
        elif k >= 3 and r % 8 == 3:
            ret = [rho(M, 3 * r), rho(M, 3 * r), 0]
        elif k == 3 and r % 8 == 5:
            ret = [2 * rho(M, r)]
        elif k >= 4 and r % 8 == 1:
            ret = [8 * rho(M, r)]
        else:
            ret = [0 for i in range(0, len(Ktil))]
        # By projection we have 
        # KW.Fg* = 1/2*(KZo - FFo).(Fg_o) = 1/2*(KZo.Fgo - Fgo^2)
        # By adjunction this is KZo.Fgo - (g - 1)
        return [
            2 * ((Ktil[i] - ret[i]) - (genera[i] - 1)) 
            for i in range(0, len(ret))
        ]


class WNr_bar:
    """The class of the surface \\bar(W)_(N,r)

    Parameters
    ----------
    N : int
        Level of the Humbert surface
    r : int
        Power of the congruence

    Attributes
    ----------
    N : int
        Level of the Humbert surface
    r : int
        Power of the congruence
    M : int
        Odd part of N
    k : int
        So that N = M * 2**k
    
    """

    def __init__(self, N, r):
        if gcd(N, r) != 1:
            raise ValueError("N,r should be coprime")
        if Ztil(N, r).kodaira_dim() == -1:
            raise ValueError(
                "W_(N,r) is rational but not defined since Z_(N,r) is rational"
            )
        self.N = N
        self.r = r
        self.k = p_adic_val(N, 2)
        self.M = N // 2 ** (self.k)

    def p_g(self):
        """BIRATIONAL INVARIANT: The geometric genus of W_(N,r).

        Returns
        -------
        int
            The geometric genus of \bar(W)_(N,r).
        
        """
        return WNr(self.N, self.r).p_g()

    def kodaira_dim(self):
        """BIRATIONAL INVARIANT: The Kodaira dimension of W_(N,r).

        Returns
        -------
        int
            The Kodaira dimension of \\bar(W)_(N,r).
        
        """
        return WNr(self.N, self.r).kodaira_dim()

    def K_Cinf(self):
        """The intersection number K.\\bar(C)_∞

        This is the intersection number of \\bar(C)_∞ with a canonical divisor
        on \\bar(W)_(N,r).

        Returns
        -------
        int
            The intersection number K.\\bar(C)_∞.

        Notes
        -----
        As per the proof of [1]_ Prop. 11.1.

        References
        ----------
        .. [1] S. Frengley, 'On the geometry of the Humbert surface of
           discriminant N^2', this article (2024).       
        
        """
        N = self.N
        r = self.r
        W = WNr(N, r)
        return W.K_Cinf() - (
            sum(
                [
                    rho(int(d), -r) * euler_phi(N // int(d))
                    for d in divisors(N)
                    if d != 1
                ]
            )
            // 2
        )


class WNr_o:
    """The class of the surface W°_(N,r)

    Parameters
    ----------
    N : int
        Level of the Humbert surface
    r : int
        Power of the congruence

    Attributes
    ----------
    N : int
        Level of the Humbert surface
    r : int
        Power of the congruence
    M : int
        Odd part of N
    k : int
        So that N = M * 2**k
    
    """

    def __init__(self, N, r):
        if gcd(N, r) != 1:
            raise ValueError("N,r should be coprime")
        if Ztil(N, r).kodaira_dim() == -1:
            raise ValueError(
                "W_(N,r) is rational but not defined since Z_(N,r) is rational"
            )
        if WNr(N, r).p_g() == 0:
            raise ValueError(
                "W^°_(N,r) is only defined for non-rational W_(N,r)"
            )
        self.N = N
        self.r = r
        self.k = p_adic_val(N, 2)
        self.M = N // 2 ** (self.k)

    def __str__(self):
        ret = "The surface W°_({N},{r})\n".format(N=self.N, r=self.r)
        ret += "κ     = {}\n".format(self.kodaira_dim())
        ret += "p_g   = {}\n".format(self.p_g())
        ret += "K^2   = {}\n".format(self.K_K())
        # ret += "χ_top = {}".format(self.chi_top())        # not implemented
        return ret

    def p_g(self):
        """BIRATIONAL INVARIANT: The geometric genus of W°_(N,r).

        Returns
        -------
        int
            The geometric genus of W°_(N,r).
        
        """
        W = WNr(self.N, self.r)
        return W.p_g()

    def kodaira_dim(self):
        """BIRATIONAL INVARIANT: The Kodaira dimension of W°_(N,r).

        Returns
        -------
        int
            The Kodaira dimension of W°_(N,r).
        
        """
        W = WNr(self.N, self.r)
        return W.kodaira_dim()

    def K_K(self):
        """Self intersection of a canonical divisor on W°_(N,r)

        Returns
        -------
        int
            The self intersection, K.K, of a canonical divisor on W°_(N,r)

        Notes
        -----
        This is just a lookup formula to [1]_ Lemma 11.2. 
        
        This is probably quite close to the K.K of the minimal model whenever
        N is coprime to 6 (I would hazard a guess that it is in all but some
        finitely many cases, maybe none at all). When N is divisible by 2 or 3,
        it's clear that we miss the mark. See [1]_ Conj. 11.7.
        
        References
        ----------
        .. [1] S. Frengley, 'On the geometry of the Humbert surface of
           discriminant N^2', this article (2024).       
        
        """
        N = self.N
        r = self.r
        M = self.M
        W = WNr(N, r)
        ret = W.K_K()
        # blow-down the components of E_(∞,d,-1)^*
        ret += (
            sum(
                [0]
                + [
                    (int(d) // 2) * rho(int(d), -r) * euler_phi(N // d)
                    for d in divisors(N)
                    if d != 1
                ]
            )
            // 2
        )
        # blow-down the components of E_(2,1)^*
        ret += s_21(N, r)
        # blow-down the components of E_(3,2)^*
        ret += s_32(N, r)
        # blow-down the exceptional F_(m,λ)*
        g0X0 = get_rat_X0_plus()
        g0X0 = [m for m in g0X0 if rho(N, r * m) != 0 if m >= 5]
        exceptional_X0 = [m for m in g0X0 if W.K_Fm(m)[0] == -1]
        ret += sum([0] + [rho(N, r * m) for m in exceptional_X0]) // 2
        # blow-down component of E_(3,1)^* which met F_(3,λ)
        ret += rho(N, 3 * r) // 2
        # and the component of E_(2,1)^* which met F_(5,λ)
        ret += rho(N, 5 * r) // 2

        # When N mod 2 == 0 need additional things
        curves_Flift_and_E2 = 0
        E31_from_Isharp = 0
        curves_Fb4 = 0
        curves_F3_b = 0
        # Above are zero except in the cases as follows
        if (N % 2 == 0):
            curves_Flift_and_E2 = rho(N // 2, 2 * r)
            curves_F3_b = rho(N, 3 * r) // 2
            if (N % 4 != 0):
                E31_from_Isharp = 2 * rho(M, r) // 4
            elif N % 4 == 0:
                if (N % 8 != 0):
                    E31_from_Isharp = rho(4, 3 * r) * rho(M, r) // 4
                    curves_Fb4 = rho(N, r) // 2
                elif (N % 8 == 0):
                    curves_Fb4 = rho(N, r) // 4
                    if N % 16 != 0:
                        E31_from_Isharp = rho(8, 5 * r) * rho(M, r) // 4
                    elif (N % 16 == 0):
                        E31_from_Isharp = rho(8, r) * rho(M, r) // 4
        ret += curves_Flift_and_E2
        ret += E31_from_Isharp
        ret += curves_Fb4
        ret += curves_F3_b

        # When N mod 9 == 3,6 need additional things
        additional = 0
        if (N % 3 == 0) and (N % 9 != 0):
            additional = rho(N // 3, r) * (4 * rho(3, r) + 5 * rho(3, 2 * r)) // 4
        ret += additional

        return ret
