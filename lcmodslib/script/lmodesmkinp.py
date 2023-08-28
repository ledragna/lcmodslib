#!/usr/bin/env python3
import os
import sys
import argparse
import json

from lcmodslib.base import gio
from lcmodslib.base import gmanip
from lcmodslib.base.utils import readlistnum

## Generic functions
#custom action argparse
def _write_json():
    """Writes a file with a template JSON file with the Gaussian Options

    """
    example = """{"qmdata": {"functional": "B3PW91",
            "basis": "gen 6D empiricaldispersion=gd3bj",
            "molchr": 0,
            "molspn": 1,
            "addroot": "freq=vcd ",
            "addline": "@SNSD.gbs",
            "mem": 10,
            "cpu": 6},
  "fragments": []
}"""

    with open('gaussian_options_example.json', 'w') as fopen:
        fopen.write(example)


class WriteTemplate(argparse.Action):
    def __init__(self, option_strings, dest, nargs=0, choices=None, const=None, **kwargs):
        super(WriteTemplate, self).__init__(option_strings=option_strings, nargs=nargs, dest=dest,**kwargs)

    def __call__(self, parser, values, namespace, option_string):
        print("   ### Printing JSON option file ###  ")
        _write_json()      
        parser.exit()


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
    parser.add_argument('--frags', action='store_true', help="add fragments number in the input (fragments must be given in the json file)")
    parser.add_argument('--silent', action='store_true', help="Suppress almost all the printing")
    parser.add_argument('--template', action=WriteTemplate,
                        help="Prints an example of the JSON with Gaussian options")
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
        optsfile = read_optfile(opts.optfile)
        qmopts = optsfile['qmdata']
    except FileNotFoundError:
        print("JSON file not found")
        sys.exit()
    natoms = len(moldata['atnum'])
    if opts.frags:
        try:
            tmpfrags = optsfile['fragments']
            if len(tmpfrags) < 2:
                print("Error: less than 2 fragments defined")
                sys.exit()
            addedindex = []
            atomfragindex = [0 for _ in range(natoms)]
            for i, frg in enumerate(tmpfrags):
                _frg = readlistnum(frg, nmax=natoms-1)
                if list(set(_frg) & set(addedindex)):
                    print("Intersections between the fragments, check the definition")
                    sys.exit()
                for x in _frg:
                    atomfragindex[x] = i
        except KeyError:
            print("No fragments in JSON")
            sys.exit()
    else:
        atomfragindex=None

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
        gio.write_gjf(hxobj.atlab, tmpgeoms, qmopts, out_file=bprefix, where=path,fragments=atomfragindex)


if __name__ == '__main__':
    main()
