import os
import sys
import argparse

from lcmodslib.base import gio
from lcmodslib.base import gmanip
# from lcmodslib.base import lmodes


# Generic functions
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
    parser.add_argument('folder', type=str,
                        help="folder where search fchk files")
    parser.add_argument('prefix', type=str,
                        help="prefix of the file names. Example of filename: prefix_bond_HO_13_12_step001_trim.fchk")
    parser.add_argument('xatoms', nargs='+', type=str,
                        help='X Atom types of XH bond')
    parser.add_argument('--quanta', type=int, default=3, choices=[1,2,3,4,5,6],
                        help="Max number of quanta")
    parser.add_argument('--nterms', type=int, default=3, choices=[1,2,3],
                        help="Max number of terms to be considered in the expansion")
    parser.add_argument('--harm', action="store_true",
                        help="Print only the Dipoles Harmonic solution")    
    parser.add_argument('--bonds', type=int, nargs='+', default=None,
                        help='extract data only for these selected bonds, a bond is specified with the two atomic indices. (e.g. --bonds 3 4 5 6 for two bonds)')
#    parser.add_argument('-w', '--where', type=str,
#                        help="whare write the input files")
    parser.add_argument('--fout', type=str, help="Output file name")
    parser.add_argument('--verboseinfo', action="store_true", help="Prints out for each bond at every step energy, apt and aat")
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
    if opts.bonds:
        if len(opts.bonds) % 2:
            print("Error. Odd number of indices")
            sys.exit()
        for atm in opts.bonds:
            if atm > len(moldata['atnum']):
                print(f'{atm} out of range')
                sys.exit()
        selbnds = [(opts.bonds[i*2]-1, opts.bonds[i*2+1]-1) for i in range(len(opts.bonds)//2)]
    else:
        selbnds = None
    # for xatm in opts.xatoms:
    hxobj = gmanip.XHstreching(moldata['atnum'],
                               moldata['atcrd'], opts.xatoms)

    if not len(hxobj.hxbonds):
        print(f"No {opts.xatoms}-H bonds found")
        sys.exit()
    if opts.prefix is None:
        prefix = "lmodes"
    else:
        prefix = opts.prefix
    print(selbnds)
    lmodesmol = gio.get_bondsdatatoobg(os.path.join(opts.folder, opts.prefix),
                                       "trim", hxobj, opts.nterms, selbnds)
    if opts.fout is None:
        fout = prefix
    else:
        fout = opts.fout

    # writes omega and chi
    with open(fout+"_omegachi.dat", "w") as fopen:
        fopen.write(lmodesmol.omgchi2string())
    hdip = False
    if opts.harm:
        hdip = True
    with open(fout+"_trns.dat", "w") as fopen:
        fopen.write(lmodesmol.lmodes2string(opts.quanta, harmdip=hdip))

    if opts.verboseinfo:
        for bnd in lmodesmol.bonds:
            tmp_data = lmodesmol.get_bondpointdata(bnd)
            # print(tmp_data['apt'])
            for key in tmp_data:
                print(f"Printing: {key} bond {bnd[0]+1}-{bnd[1]+1}")
                fname = f"{fout}_bond_{bnd[0]+1}_{bnd[1]+1}_{key}_dz.dat"
                with open(fname, "w") as fopen:
                    fopen.write(gio.write_dataarray(tmp_data[key]))
            

if __name__ == '__main__':
    main()
