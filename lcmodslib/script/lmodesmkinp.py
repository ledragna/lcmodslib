#!/usr/bin/env python3
"""
Example of Json
{"qmdata": {"functional": "B3PW91",
                    "basis": "gen 6D 10F empiricaldispersion=gd3bj",
                    "molchr": 0,
                    "molspn": 1,
                    "addroot": "freq=vcd ",
                    "addline": "@/home/m.fuse/basis/SNSD.gbs"},
                    "mem": 20,
                    "cpu": 12
                    }
"""
import os
import sys
import argparse
import copy
import json
import warnings

from estampes.data.atom import atomic_data

from lcmodslib.base import gio
from lcmodslib.base import gmanip

## Generic functions
def read_optfile(fname):
    """
    Read the option file. expected in json style
    """
    if not os.path.exists(fname):
        raise FileNotFoundError('Missing JSON file')
    with open(fname, 'r', encoding='UTF-8') as fopen:
        data = json.load(fopen)
    return data

def build_parser():
    """Builds options parser.
    Builds the full option parser.
    Returns
    -------
    :obj:`ArgumentParser`
        `ArgumentParser` object
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('fname', type=str, help="The gaussian log file with the equilibrium geometry")
    parser.add_argument('optfile', nargs='?', help='option file in json with QM options')
    parser.add_argument('xatoms', nargs='+', type=str,
                        help='X Atom types of XH bond')
    parser.add_argument('--minpos', type=float, default=-0.33,
                        help="Minimum position in angstrom")
    parser.add_argument('--nstep', type=int, default=13,
                        help="Number of steps from the minimum position")
    parser.add_argument('--stepsize', type=float, default=0.064,
                        help="step size in angstrom")
    parser.add_argument('-w', '--where', type=str,
                        help="whare write the input files")
    parser.add_argument('--prefix', type=str,
                        help="prefix to add to the input files")
    parser.add_argument('--silent', action='store_true', help="Suppress almost all the printing")
    return parser


def main():
    par = build_parser()
    opts = par.parse_args()
    try:
        moldata = gio.get_mol_data(opts.fname)
    except FileNotFoundError:
        print(f"{opts.fname} Not Found")
        sys.exit()
    try:
        qmopts = read_optfile(opts.optfile)['qmdata']
    except FileNotFoundError:
        print("JSON file not found")
        sys.exit()
    if not opts.where:
        path = "."
    else:
        path = opts.where
    if not os.path.exists(path):
        print(f'{path} not exist. Created.')
        os.mkdir(path)

    hxobj = gmanip.XHstreching(moldata['atnum'],
                               moldata['atcrd'], opts.xatoms)
    if not len(hxobj.hxbonds):
        print(f"No {opts.xatoms}-H bonds found")
        sys.exit()
    if opts.prefix is None:
        prefix = "lmodes"
    else:
        prefix = opts.prefix

    lbonds = hxobj.hxbonds
    xtps = hxobj.getsecatom()
    for i, bond in enumerate(lbonds):
        bprefix = prefix + f"_bond_H{xtps[i]:s}_{bond[0]+1:02d}_{bond[1]+1:02d}"
        tmpgeoms = hxobj.scanhx(bond, lower=opts.minpos,
                                step=opts.stepsize, nstep=opts.nstep)
        gio.write_gjf(hxobj.atlab, tmpgeoms, qmopts, out_file=bprefix, where=path)




if __name__ == '__main__':
    main()