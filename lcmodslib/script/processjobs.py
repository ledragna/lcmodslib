import os
import sys
import argparse
import copy
import json
import warnings

from estampes.data.atom import atomic_data

from lcmodslib.base import gio
from lcmodslib.base import gmanip
from lcmodslib.base import lmodes

## Generic functions
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
    parser.add_argument('xatoms', nargs='?', default='C', type=str,
                        help='X Atom type of XH bond')
    parser.add_argument('folder', type=str,
                        help="folder where search fchk files")
    parser.add_argument('prefix', type=str,
                        help="prefix of the file names")
    parser.add_argument('--quanta', type=int, default=3, choices=[1,2,3],
                        help="Max number of quanta")
    parser.add_argument('--nterms', type=int, default=3, choices=[1,2,3],
                        help="Max number of terms to be considered in the expansion")
    parser.add_argument('-w', '--where', type=str,
                        help="whare write the input files")
    parser.add_argument('--prefix', type=str,
                        help="prefix to add to the input files")
    parser.add_argument('--fout', type=str, help="Output file name")
    return parser


def main():
    par = build_parser()
    opts = par.parse_args()
    try:
        moldata = gio.get_mol_data(opts.fname)
    except FileNotFoundError:
        print(f"{opts.fname} Not Found")
        sys.exit()
    if not os.path.exists(opts.folder):
        print(f'{opts.folder} not exist.')
        sys.exit()
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
    lmodesmol = gio.get_bondsdatatoobg(os.path.join(opts.folder, opts.prefix),
                                       opts.xatoms, "trim", hxobj.atnum,
                                       hxobj.refcrd, opts.nterms, *lbonds)
    if opts.fout is None:
        fout = prefix
    else:
        fout = self.fout

    # writes omega and chi
    with open(fout+"_omegachi.dat", "w") as fopen:
        fopen.write(lmodesmol.omgchi2string())
    with open(fout+"_trns.dat", "w") as fopen:
        fopen.write(lmodesmol.lmodes2string(opts.quanta))

if __name__ == '__main__':
    main()