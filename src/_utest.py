#!/usr/bin/env python3
import os
import re
import argparse
import subprocess

DEPEND = {
    "basis.f": ["eval.f", ],
    "rhf.f": ["mtrx.f", "eval.f", "basis.f", ],
    "uhf.f": ["mtrx.f", "eval.f", "basis.f", ],
    "rohf.f": ["mtrx.f", "eval.f", "basis.f", "uhf.f", ],
    "ci.f": ["mtrx.f", "eval.f", "basis.f", "uga.f", ],
    }
EXE = "unitest.x"
CMD = "{fcc} {fcflags} {depends} {target} _test/test_{target} -o {exe}"


def main():
    """Unit test driver routine. """
    args = getargs()
    keys = {
        "fcc": args.fcc,
        "fcflags": args.fcflags,
        "depends": "",
        "target": args.target,
        "exe": EXE,
        }
    # Resolve dependence.
    if args.target in DEPEND:
        keys["depends"] = " ".join(DEPEND[args.target])
    cmd = re.sub(r"\s+", " ", CMD.format(**keys))
    print(cmd)
    try:
        subprocess.check_call(cmd.split())
    except subprocess.CalledProcessError:
        exit(1)
    try:
        cmd = './{0}'.format(EXE)
        subprocess.call(cmd.split())
    finally:
        os.remove(EXE)


def getargs():
    """ """
    fcflags=" ".join([
        "-O2",
        "-lblas", "-llapack",
        "-g", "-fbounds-check",
        "-fopenmp",
        "-fdefault-integer-8",
        # "-Wall", 
        ])
    parser = argparse.ArgumentParser(
        prog="_unittest",
        description="Execute unit tests.",
        )
    parser.add_argument(
        "target",
        help="target *.f source code",
        )
    parser.add_argument(
        "--fcc",
        help="Fortran compiler (default: 'gfortran')",
        default='gfortran',
        )
    parser.add_argument(
        "--fcflags",
        help="Fortran compiler (default: '{0}')".format(fcflags),
        default=fcflags,
        )
    return parser.parse_args()


if __name__ == "__main__":
    main()
