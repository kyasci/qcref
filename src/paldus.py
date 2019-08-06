#!/usr/bin/env python3
import numpy as np

__all__ = ["Paldus", ]


class Paldus0(object):
    """ """
    def __init__(self, nea, nac, mult, nx):
        """ """
        # Paldus table for (CAS).
        self._tabcas = paldus(nea, nac, mult)
        # Operation list (CAS).
        self._config = config(self._tabcas)
        # CAS <-> RAS index.
        mcr, mrc = cr_and_rc(
        cf = np.array(self._config)
        mcr = []
        mrc = [0 for i in self._tabcas]
        for i, cfi in enumerate(cf):
            ndiff = np.sum(np.abs(cfi-cf[0]))//2
            if ndiff <= nx:
                mrc[i] = len(mcr)
                mcr.append(i)
        print(mcr)
        print(mrc)

    def get_tabcas(self, ):
        """ """
        return self._tabcas

    def get_confcas(self, ):
        """ """
        return self._config


def main():
    """ """
    nea = 3
    nac = 3
    mult = 2
    nx = 2
    paldus = Paldus(nea, nac, mult, nx)
    # tabcas = paldus.get_tabcas()
    # config = np.array(paldus.get_confcas())
    # tab = paldus.get_tab()
    # print(len(tabcas))
    # print(len(tab))


def paldus(nea, nac, mult):
    """ """
    mb = mult - 1
    ma = (nea - mb) // 2
    mc = nac - ma - mb
    ls = [[[ma, mb, mc]]]
    km = np.array([[0, 0, 1], [0, 1, 0], [1, -1, 1], [1, 0, 0], ])
    lpaldus = [[[ma, mb, mc]]]
    for i in range(nac - 1):
        ltmp = []
        for l in lpaldus:
            for k in km:
                if proceed(l, k):
                    lt = l[:]
                    lt.append((np.array(l[-1])-np.array(k)).tolist())
                    ltmp.append(lt)
        lpaldus[:] = ltmp[:]
    return lpaldus


def config(lpaldus):
    """ """
    lop = []
    for l in lpaldus:
        s = []
        for i, x in enumerate(l[:-1]):
            diff = (np.array(l[i]) - np.array(l[i+1])).tolist()
            s.append(diff_op(*diff))
        s.append(diff_op(*l[-1]))
        lop.append(s[::-1])
    return lop


def diff_op(ka, kb, kc):
    """ """
    if ka == 0 and kb == 0 and kc == 1:
        ms = 0
    if ka == 0 and kb == 1 and kc == 0:
        ms = 1
    if ka == 1 and kb == -1 and kc == 1:
        ms = -1
    if ka == 1 and kb == 0 and kc == 0:
        ms = 2
    return ms


def proceed(l, k):
    """ """
    res = True
    for li, ki in zip(l[-1], k):
        if li - ki < 0:
            res = False
            break
    return res


if __name__ == "__main__":
    main()
