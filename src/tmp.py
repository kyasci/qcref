#!/usr/bin/env python3
import numpy as np


class Paldus(object):
    """ """
    def __init__(self, nea, nac, mult, nx):
        """ """
        self._nx = nx
        self._paltab = _eval_paltab(nea, nac, mult)

    def get_tabcas(self, ):
        """ """
        return self._paltab

    def get_config(self, ):
        """ """
        config = []
        for l in self._paltab:
            s = []
            for i, x in enumerate(l[:-1]):
                diff = (np.array(l[i]) - np.array(l[i+1])).tolist()
                s.append(_conf_from(diff))
            s.append(_conf_from(l[-1]))
            config.append(s[::-1])
        return config

    def get_mrc(self, ):
        """RAS -> CAS index."""
        cf = np.array(self.get_config())
        mrc = []
        for i, cfi in enumerate(cf):
            ndiff = np.sum(np.abs(cfi-cf[0]))//2
            if ndiff <= self._nx:
                mrc.append(i)
        return mrc

    def get_mcr(self, ):
        """CAS -> RAS index."""
        cf = np.array(self.get_config())
        mcr = [-1 for i in self._paltab]
        npt = 0
        for i, cfi in enumerate(cf):
            ndiff = np.sum(np.abs(cfi-cf[0]))//2
            if ndiff <= self._nx:
                mcr[i] = npt
                npt += 1
        return mcr



def _conf_from(diff):
    """ """
    ka, kb, kc = diff
    if ka == 0 and kb == 0 and kc == 1:
        ms = 0
    if ka == 0 and kb == 1 and kc == 0:
        ms = 1
    if ka == 1 and kb == -1 and kc == 1:
        ms = -1
    if ka == 1 and kb == 0 and kc == 0:
        ms = 2
    return ms


def _eval_paltab(nea, nac, mult):
    """Evaluate the Paldus table."""
    # The first row.
    mb = mult - 1
    ma = (nea - mb) // 2
    mc = nac - ma - mb
    ls = [[[ma, mb, mc]]]
    # Four possible operations.
    km = np.array([[0, 0, 1], [0, 1, 0], [1, -1, 1], [1, 0, 0], ])
    # Evaluate the Paldus table.
    lpaldus = [[[ma, mb, mc]]]
    for i in range(nac - 1):
        ltmp = []
        for l in lpaldus:
            for k in km:
                if all(np.array(l[-1]) - np.array(k) >= 0):
                    lt = l[:]
                    lt.append((np.array(l[-1])-np.array(k)).tolist())
                    ltmp.append(lt)
        lpaldus[:] = ltmp[:]
    return lpaldus


def main():
    """ """
    nea, nac, mult = 3, 3, 2
    # nea, nac, mult = 4, 4, 1
    nx = 1
    paldus = Paldus(nea, nac, mult, nx)
    paltab = paldus.get_tabcas()
    conf = paldus.get_config()
    mcr = paldus.get_mcr()
    mrc = paldus.get_mrc()
    for i, pt in enumerate(paltab):
        if mcr[i] >= 0:
            # print(pt)
            print(conf[i])


if __name__ == "__main__":
    main()
