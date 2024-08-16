"""
This is a script which generates in ./figs/N-r/ pictures in a neighbourhood of
the cusps on all Z_{N,r} for which W_{N,r} is either a K3 or properly elliptic
surface.
"""

import matplotlib.pyplot as plt
from humbert.znr import Ztil
from humbert.znr import WNr


def main():
    """
    Main
    """
    for N in range(15, 25):
        for r in possible_r(N):
            W = WNr(N, r)
            if W.p_g() in [1, 2]:
                Z = Ztil(N, r)
                cc = Z.sketch_cusps(
                    tex=True, compile_tex=True, open_pdf=True, disp=False
                )
                plt.close("all")


if __name__ == "__main__":
    main()
