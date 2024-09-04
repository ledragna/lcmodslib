import os
import numpy as np
import glob
import typing as tp
from estampes.parser.base import DataFile
from estampes.base import QLabel
from estampes.tools.atom import convert_labsymb
from estampes.data.physics import PHYSFACT
from estampes.data.atom import atomic_data
from lcmodslib.base.lmodes import LocalModes, LmodDeriv


def get_mol_data(fname):
    dkeys = {'Energy': QLabel(quantity=1),
             'atcrd': QLabel(quantity='atcrd', descriptor='last'),
             'atnum': QLabel(quantity='atnum'),
             }

    dfile = DataFile(fname)
    res = {}
    res['fname'] = fname
    data = dfile.get_data(**dkeys)
    res['atnum'] = data['atnum'].data
    res['atlab'] = convert_labsymb(True, *data['atnum'].data)
    res['atcrd'] = np.array(data['atcrd'].data)*PHYSFACT.bohr2ang
    res['eng'] = data['Energy'].data
    atdat = atomic_data(*set(res['atlab']))
    res['atmass'] = np.array([atdat[at]['mass'] for at in res['atlab']])
    return res


def write_xyz(atnum, comments, *geoms):
    """Writes molecular geometries into an traj xyz file

    Args:
        atnum (list(int)): Atomic numbers
        comments (list(str)): List of comments to be added to the geometries
        geoms (np.array): the Cartesian coordinates to be written

    Raises:
        IndexError: _description_

    Returns:
        str: the string to be written in a file
    """
    line = '{at:4s}{xyz[0]:12.6f}{xyz[1]:12.6f}{xyz[2]:12.6f}\n'
    # ngeom = len(geoms)
    natoms = len(atnum)
    string = ""
    if isinstance(comments, (list, tuple, np.ndarray)):
        ref = comments[0]
        labels = comments
    else:
        ref = comments
        labels = [comments]
    if len(geoms) == 1:
        fmt = '{val}\n'
    else:
        fmt = 'geom:{gid:03d} comm: {val}\n'
    if len(labels) != len(geoms):
        if len(labels) == 1:
            labels = [ref for _ in range(len(geoms))]
        else:
            raise IndexError('Mismatch between geometries and labels')
    for i, geom in enumerate(geoms):
        string += '{:d}\n'.format(natoms)
        string += fmt.format(gid=i+1, val=labels[i])
        for iat, xyz in enumerate(geom):
            string += line.format(at=atnum[iat], xyz=xyz)
    return string.strip()


def write_gjf(atlab, geoms, qmdata, out_file='sel_frame',
              where="", fragments=None):
    """Write a gaussian input file
    Args:
        atlabs ([type]): [description]
        geoms ([type]): [description]
        out_file (str, optional): [description]. Defaults to 'sel_frame'.
        where (str): where save the files
        fragments (list(int)): list of fragment indices (default=None)
    """
    template = """%mem={MEM}GB
%nprocshared={CPU}
%chk={CHK}
#p {FUN}/{BASIS}
 nosymm {ADDROOT}

{COMM}

{MOLCHRSPN}
"""
    line = '{at:>4s}{add:20s}{xyz[0]:12.6f}{xyz[1]:12.6f}{xyz[2]:12.6f}\n'
    # for each structure
    tmplcomm = '{val}'
    if len(geoms) > 1:
        tmplcomm = 'geom:{gid:03d} {val}'

    for i, geom in enumerate(geoms):
        geomstr = ""
        chgspn = ""
        if fragments:
            if (len(qmdata['molchr']) != len(qmdata['molspn'])) or (len(qmdata['molchr']) != len(list(set(fragments)))+1):
                raise TypeError("Error in fragments definition")
            for atindx in range(len(qmdata['molchr'])):
                chgspn += f"{qmdata['molchr'][atindx]} {qmdata['molspn'][atindx]} "
        else:
            chgspn = f"{qmdata['molchr']} {qmdata['molspn']}"
        for iat, xyz in enumerate(geom):
            if fragments:
                add = f"(Fragment={fragments[iat]+1})"
                # if (len(qmdata['molchr']) != len(qmdata['molspn'])) or (len(qmdata['molchr']) != len(list(set(fragments)))+1):
                #    raise TypeError("Error in fragments definition")
                # chgspn = ""
                for atindx in range(len(qmdata['molchr'])):
                    chgspn += f"{qmdata['molchr'][atindx]} {qmdata['molspn'][atindx]} "
            else:
                add = ""
                # chgspn = f"{qmdata['molchr']} {qmdata['molspn']}"
            geomstr += line.format(at=atlab[iat], add=add, xyz=xyz)
        fout = out_file+'_step{:03d}.gjf'.format(i)
        actout = os.path.join(where, fout)
        with open(actout, 'w') as fopen:
            fopen.write(template.format(MEM=qmdata['mem'],
                                        CPU=qmdata['cpu'],
                                        CHK=fout[:-3]+'chk',
                                        FUN=qmdata['functional'],
                                        BASIS=qmdata['basis'],
                                        MOLCHRSPN=chgspn,
                                        ADDROOT=qmdata['addroot'],
                                        COMM=tmplcomm.format(gid=i+1, val="test")))
            fopen.write(geomstr)
            fopen.write('\n')
            fopen.write(f'{qmdata["addline"]}\n\n')


def get_fchk(fname):
    dkeys2 = {'eng': QLabel(quantity=1),
              'atcrd': QLabel(quantity='atcrd', descriptor='last'),
              'atnum': QLabel(quantity='atnum'),
              }
    # print(fname)
    dfile = DataFile(fname)
    res = {}
    res['fname'] = fname
    data = dfile.get_data(**dkeys2)
    atmnum = len(data['atnum'].data)
    res['atnum'] = data['atnum'].data
    res['atlab'] = convert_labsymb(True, *data['atnum'].data)
    res['atcrd'] = np.array(data['atcrd'].data)*PHYSFACT.bohr2ang
    res['eng'] = data['eng'].data

    try:
        dkey = {'apt': QLabel(quantity=101, derorder=1)}
        data = dfile.get_data(**dkey)
        res['apt'] = np.array(data['apt'].data).reshape(atmnum, 3, 3)
    except:
        res['apt'] = None
        # res['aat'] = None
    try:
        dkey = {'aat': QLabel(quantity=102, derorder=1)}
        data = dfile.get_data(**dkey)
        res['aat'] = np.array(data['aat'].data).reshape(atmnum, 3, 3)
    except:
        res['aat'] = None

    return res


def get_bondsdatatoobg(prefix, suffix, hxobj, nterms, selbnds=None):
    res = LocalModes(hxobj.atnum, hxobj.refcrd, nterms)
    sys_bonds = hxobj.hxbonds
    print(sys_bonds)
    if not selbnds:
        bonds = sys_bonds
    else:
        bonds = []
        for sbnd in selbnds:
            if sbnd in sys_bonds:
                bonds.append(sbnd)
            elif (sbnd[1], sbnd[0]) in sys_bonds:
                bonds.append((sbnd[1], sbnd[0]))
            else:
                print(f"{sbnd[0]+1}-{sbnd[1]+1} is not a bond of the selected types. Skipped" )
        if not bonds:
            raise BaseException("No bonds")

    fname = prefix +"_bond_H{}_{:02d}_{:02d}_step"
    fname2 = prefix +"_bond_H{}_{:02d}_{:02d}_step{:03d}_"+suffix+".fchk"
    for bnd in bonds:
        # print(bnd)
        atype = hxobj.getsecatom(bnd)
        lfiles = glob.glob(fname.format(atype, bnd[0]+1, bnd[1]+1)+"*.fchk")
        # print(lfiles)
        tmpres = {'eng': [], 'len': [],
                  'apt1': [], 'aat1': [],
                  'apt2': [], 'aat2': []}
        for i in range(len(lfiles)):
            # print(i)
            tmp_data = get_fchk(fname2.format(atype, bnd[0]+1, bnd[1]+1, i))
            # print(tmp_data)
            tmpres['eng'].append(tmp_data['eng'])
            tmpres['len'].append(np.abs(tmp_data['atcrd'][bnd[0], 2]-tmp_data['atcrd'][bnd[1], 2]))
            if not tmp_data['apt'] is None:
                tmpres['apt1'].append(tmp_data['apt'][bnd[0], 2, :])
                tmpres['apt2'].append(tmp_data['apt'][bnd[1], 2, :])
            else:
                tmpres['apt1'].append(None)
                tmpres['apt2'].append(None)
            if not tmp_data['aat'] is None:
                tmpres['aat1'].append(tmp_data['aat'][bnd[0], 2, :])
                tmpres['aat2'].append(tmp_data['aat'][bnd[1], 2, :])
            else:
                tmpres['aat1'].append(None)
                tmpres['aat2'].append(None)

        tmp_bond = LmodDeriv()
        tmp_bond.bond = bnd
        tmp_bond.eng = tmpres['eng']
        tmp_bond.blen = tmpres['len']
        tmp_bond.apt = (tmpres['apt1'], tmpres['apt2'])
        tmp_bond.aat = (tmpres['aat1'], tmpres['aat2'])
        res.addbond(tmp_bond)
    return res


def write_dataarray(data: tp.Tuple[np.ndarray,
                                   tp.Union[np.ndarray, list]]) -> str:
    """Writes a data array in a string format

    Args:
        data (tp.Tuple[np.ndarray, np.ndarray]): The data array to be written

    Returns:
        str: The string to be written
    """
    if isinstance(data[1], list):
        _tmpdata = np.hstack((data[1][0], data[1][1]))
    else:
        _tmpdata = data[1]
    if data[0].shape[0] != _tmpdata.shape[0]:
        raise ValueError("The two arrays have different sizes")
    string = ""
    if _tmpdata.ndim > 1:
        nval = _tmpdata[0].flatten().shape[0]
    else:
        nval = 1
    tmpl = "{:12.6f} " + "{:18.10E} " * (nval) + "\n"
    for i in range(data[0].shape[0]):
        string += tmpl.format(data[0][i], *_tmpdata[i].flatten())
    return string
