"""
classnumbers.py
===============

This module consists of a single global variable which hardcodes the class
numbers of imaginary quadratic orders of discriminant at most -1000.

Functions
---------
- class_number(d)

"""

CLASS_NUMBER_DICT = {
    -3: 1,
    -4: 1,
    -7: 1,
    -8: 1,
    -11: 1,
    -12: 1,
    -15: 2,
    -16: 1,
    -19: 1,
    -20: 2,
    -23: 3,
    -24: 2,
    -27: 1,
    -28: 1,
    -31: 3,
    -32: 2,
    -35: 2,
    -36: 2,
    -39: 4,
    -40: 2,
    -43: 1,
    -44: 3,
    -47: 5,
    -48: 2,
    -51: 2,
    -52: 2,
    -55: 4,
    -56: 4,
    -59: 3,
    -60: 2,
    -63: 4,
    -64: 2,
    -67: 1,
    -68: 4,
    -71: 7,
    -72: 2,
    -75: 2,
    -76: 3,
    -79: 5,
    -80: 4,
    -83: 3,
    -84: 4,
    -87: 6,
    -88: 2,
    -91: 2,
    -92: 3,
    -95: 8,
    -96: 4,
    -99: 2,
    -100: 2,
    -103: 5,
    -104: 6,
    -107: 3,
    -108: 3,
    -111: 8,
    -112: 2,
    -115: 2,
    -116: 6,
    -119: 10,
    -120: 4,
    -123: 2,
    -124: 3,
    -127: 5,
    -128: 4,
    -131: 5,
    -132: 4,
    -135: 6,
    -136: 4,
    -139: 3,
    -140: 6,
    -143: 10,
    -144: 4,
    -147: 2,
    -148: 2,
    -151: 7,
    -152: 6,
    -155: 4,
    -156: 4,
    -159: 10,
    -160: 4,
    -163: 1,
    -164: 8,
    -167: 11,
    -168: 4,
    -171: 4,
    -172: 3,
    -175: 6,
    -176: 6,
    -179: 5,
    -180: 4,
    -183: 8,
    -184: 4,
    -187: 2,
    -188: 5,
    -191: 13,
    -192: 4,
    -195: 4,
    -196: 4,
    -199: 9,
    -200: 6,
    -203: 4,
    -204: 6,
    -207: 6,
    -208: 4,
    -211: 3,
    -212: 6,
    -215: 14,
    -216: 6,
    -219: 4,
    -220: 4,
    -223: 7,
    -224: 8,
    -227: 5,
    -228: 4,
    -231: 12,
    -232: 2,
    -235: 2,
    -236: 9,
    -239: 15,
    -240: 4,
    -243: 3,
    -244: 6,
    -247: 6,
    -248: 8,
    -251: 7,
    -252: 4,
    -255: 12,
    -256: 4,
    -259: 4,
    -260: 8,
    -263: 13,
    -264: 8,
    -267: 2,
    -268: 3,
    -271: 11,
    -272: 8,
    -275: 4,
    -276: 8,
    -279: 12,
    -280: 4,
    -283: 3,
    -284: 7,
    -287: 14,
    -288: 4,
    -291: 4,
    -292: 4,
    -295: 8,
    -296: 10,
    -299: 8,
    -300: 6,
    -303: 10,
    -304: 6,
    -307: 3,
    -308: 8,
    -311: 19,
    -312: 4,
    -315: 4,
    -316: 5,
    -319: 10,
    -320: 8,
    -323: 4,
    -324: 6,
    -327: 12,
    -328: 4,
    -331: 3,
    -332: 9,
    -335: 18,
    -336: 8,
    -339: 6,
    -340: 4,
    -343: 7,
    -344: 10,
    -347: 5,
    -348: 6,
    -351: 12,
    -352: 4,
    -355: 4,
    -356: 12,
    -359: 19,
    -360: 8,
    -363: 4,
    -364: 6,
    -367: 9,
    -368: 6,
    -371: 8,
    -372: 4,
    -375: 10,
    -376: 8,
    -379: 3,
    -380: 8,
    -383: 17,
    -384: 8,
    -387: 4,
    -388: 4,
    -391: 14,
    -392: 8,
    -395: 8,
    -396: 6,
    -399: 16,
    -400: 4,
    -403: 2,
    -404: 14,
    -407: 16,
    -408: 4,
    -411: 6,
    -412: 5,
    -415: 10,
    -416: 12,
    -419: 9,
    -420: 8,
    -423: 10,
    -424: 6,
    -427: 2,
    -428: 9,
    -431: 21,
    -432: 6,
    -435: 4,
    -436: 6,
    -439: 15,
    -440: 12,
    -443: 5,
    -444: 8,
    -447: 14,
    -448: 4,
    -451: 6,
    -452: 8,
    -455: 20,
    -456: 8,
    -459: 6,
    -460: 6,
    -463: 7,
    -464: 12,
    -467: 7,
    -468: 8,
    -471: 16,
    -472: 6,
    -475: 4,
    -476: 10,
    -479: 25,
    -480: 8,
    -483: 4,
    -484: 6,
    -487: 7,
    -488: 10,
    -491: 9,
    -492: 6,
    -495: 16,
    -496: 6,
    -499: 3,
    -500: 10,
    -503: 21,
    -504: 8,
    -507: 4,
    -508: 5,
    -511: 14,
    -512: 8,
    -515: 6,
    -516: 12,
    -519: 18,
    -520: 4,
    -523: 5,
    -524: 15,
    -527: 18,
    -528: 8,
    -531: 6,
    -532: 4,
    -535: 14,
    -536: 14,
    -539: 8,
    -540: 6,
    -543: 12,
    -544: 8,
    -547: 3,
    -548: 8,
    -551: 26,
    -552: 8,
    -555: 4,
    -556: 9,
    -559: 16,
    -560: 12,
    -563: 9,
    -564: 8,
    -567: 12,
    -568: 4,
    -571: 5,
    -572: 10,
    -575: 18,
    -576: 8,
    -579: 8,
    -580: 8,
    -583: 8,
    -584: 16,
    -587: 7,
    -588: 6,
    -591: 22,
    -592: 4,
    -595: 4,
    -596: 14,
    -599: 25,
    -600: 8,
    -603: 4,
    -604: 7,
    -607: 13,
    -608: 12,
    -611: 10,
    -612: 8,
    -615: 20,
    -616: 8,
    -619: 5,
    -620: 12,
    -623: 22,
    -624: 8,
    -627: 4,
    -628: 6,
    -631: 13,
    -632: 8,
    -635: 10,
    -636: 10,
    -639: 14,
    -640: 8,
    -643: 3,
    -644: 16,
    -647: 23,
    -648: 6,
    -651: 8,
    -652: 3,
    -655: 12,
    -656: 16,
    -659: 11,
    -660: 8,
    -663: 16,
    -664: 10,
    -667: 4,
    -668: 11,
    -671: 30,
    -672: 8,
    -675: 6,
    -676: 6,
    -679: 18,
    -680: 12,
    -683: 5,
    -684: 12,
    -687: 12,
    -688: 6,
    -691: 5,
    -692: 14,
    -695: 24,
    -696: 12,
    -699: 10,
    -700: 6,
    -703: 14,
    -704: 12,
    -707: 6,
    -708: 4,
    -711: 20,
    -712: 8,
    -715: 4,
    -716: 15,
    -719: 31,
    -720: 8,
    -723: 4,
    -724: 10,
    -727: 13,
    -728: 12,
    -731: 12,
    -732: 8,
    -735: 16,
    -736: 8,
    -739: 5,
    -740: 16,
    -743: 21,
    -744: 12,
    -747: 6,
    -748: 6,
    -751: 15,
    -752: 10,
    -755: 12,
    -756: 12,
    -759: 24,
    -760: 4,
    -763: 4,
    -764: 13,
    -767: 22,
    -768: 8,
    -771: 6,
    -772: 4,
    -775: 12,
    -776: 20,
    -779: 10,
    -780: 12,
    -783: 18,
    -784: 8,
    -787: 5,
    -788: 10,
    -791: 32,
    -792: 8,
    -795: 4,
    -796: 9,
    -799: 16,
    -800: 12,
    -803: 10,
    -804: 12,
    -807: 14,
    -808: 6,
    -811: 7,
    -812: 12,
    -815: 30,
    -816: 12,
    -819: 8,
    -820: 8,
    -823: 9,
    -824: 20,
    -827: 7,
    -828: 6,
    -831: 28,
    -832: 8,
    -835: 6,
    -836: 20,
    -839: 33,
    -840: 8,
    -843: 6,
    -844: 9,
    -847: 10,
    -848: 12,
    -851: 10,
    -852: 8,
    -855: 16,
    -856: 6,
    -859: 7,
    -860: 14,
    -863: 21,
    -864: 12,
    -867: 6,
    -868: 8,
    -871: 22,
    -872: 10,
    -875: 10,
    -876: 12,
    -879: 22,
    -880: 8,
    -883: 3,
    -884: 16,
    -887: 29,
    -888: 12,
    -891: 6,
    -892: 7,
    -895: 16,
    -896: 16,
    -899: 14,
    -900: 8,
    -903: 16,
    -904: 8,
    -907: 3,
    -908: 15,
    -911: 31,
    -912: 8,
    -915: 8,
    -916: 10,
    -919: 19,
    -920: 20,
    -923: 10,
    -924: 12,
    -927: 20,
    -928: 4,
    -931: 6,
    -932: 12,
    -935: 28,
    -936: 12,
    -939: 8,
    -940: 6,
    -943: 16,
    -944: 18,
    -947: 5,
    -948: 12,
    -951: 26,
    -952: 8,
    -955: 4,
    -956: 15,
    -959: 36,
    -960: 8,
    -963: 6,
    -964: 12,
    -967: 11,
    -968: 10,
    -971: 15,
    -972: 9,
    -975: 16,
    -976: 12,
    -979: 8,
    -980: 12,
    -983: 27,
    -984: 12,
    -987: 8,
    -988: 6,
    -991: 17,
    -992: 16,
    -995: 8,
    -996: 12,
    -999: 24,
}


def class_number(d):
    """
    The class number of the imaginary quadratic order of discriminant 
    -1000 < d < 0.

    Parameters
    ----------
    d : int
        A discriminant `d` < 0, that is an integer 0 or 1 mod 4, and has 
        absolute value at most 1000.

    Returns
    -------
    int
        Returns the class number of the quadratic order of discriminant `d`. 
        Equivalently, the number of primitive binary quadratic forms of 
        discriminant `d`.

    Notes
    -----
    This is not an implementation of Gauss reduction, it is a lookup table for 
    small discriminants.
    """
    if d > 0:
        raise ValueError(
            "Only implemented imaginary quadratic fields (i.e., d < 0)"
        )
    elif d <= -1000:
        raise NotImplementedError(
            "Only implemented for -1000 < d < 0"
        )
    return CLASS_NUMBER_DICT[d]
