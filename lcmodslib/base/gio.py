import os
import numpy as np
import glob
from estampes.parser import DataFile, build_qlabel
from estampes.tools.atom import convert_labsymb
from estampes.data.physics import PHYSFACT
from estampes.data.atom import atomic_data
from lcmodslib.base.lmodes import LocalModes, LmodDeriv

def get_mol_data(fname):
    dkeys = {'Energy': build_qlabel(1),
             'atcrd': build_qlabel('atcrd', 'last'),
             'atnum': build_qlabel('atnum'),
            }

    dfile = DataFile(fname)
    res = {}
    res['fname'] = fname
    data = dfile.get_data(*dkeys.values())
    res['atnum'] = data[dkeys['atnum']]['data']
    res['atlab'] = convert_labsymb(True, *data[dkeys['atnum']]['data'])
    res['atcrd'] = np.array(data[dkeys['atcrd']]['data'])*PHYSFACT.bohr2ang
    res['eng'] = data[dkeys['Energy']]['data']
    atdat = atomic_data(*set(res['atlab']))
    res['atmass'] = np.array([atdat[at]['mass'] for at in res['atlab']])
    return res


def write_xyz(atnum, comments, *geoms):
    line = '{at:4s}{xyz[0]:12.6f}{xyz[1]:12.6f}{xyz[2]:12.6f}\n'
    ngeom = len(geoms)
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

def write_gjf(atlab, geoms, qmdata, out_file='sel_frame', where=""):
    """Write a gaussian input file
    Args:
        atlabs ([type]): [description]
        geoms ([type]): [description]
        out_file (str, optional): [description]. Defaults to 'sel_frame'.
    """
    template ="""%mem={MEM}GB
%nprocshared={CPU}
%chk={CHK}
#p {FUN}/{BASIS}
 nosymm {ADDROOT}
 
{COMM}

{MOLCHR} {MOLSPN}
"""
    line = '{at:4s}{xyz[0]:12.6f}{xyz[1]:12.6f}{xyz[2]:12.6f}\n'
    # for each structure
    tmplcomm = '{val}'
    if len(geoms) > 1:
        tmplcomm = 'geom:{gid:03d} {val}'
        
    
    for i, geom in enumerate(geoms):
        geomstr = ""
        for iat, xyz in enumerate(geom):
            geomstr += line.format(at=atlab[iat], xyz=xyz)
        fout = out_file+'_step{:03d}.gjf'.format(i)
        actout = os.path.join(where,fout)
        with open(actout, 'w') as fopen:
            fopen.write(template.format(MEM=qmdata['mem'],
                                        CPU=qmdata['cpu'],
                                        CHK=fout[:-3]+'chk',
                                        FUN=qmdata['functional'],
                                        BASIS=qmdata['basis'],
                                        MOLCHR=qmdata['molchr'],
                                        MOLSPN=qmdata['molspn'],
                                        ADDROOT=qmdata['addroot'],
                                        COMM=tmplcomm.format(gid=i+1, val="test")))
            fopen.write(geomstr)
            fopen.write('\n')
            fopen.write(f'{qmdata["addline"]}\n\n')
        
def get_fchk(fname):
    dkeys = {'eng': build_qlabel(1),
             'atcrd': build_qlabel('atcrd', 'last'),
             'atnum': build_qlabel('atnum'),
             'aat': build_qlabel(102,None,1),
             'apt': build_qlabel(101,None,1)
            }
    dkeys2 = {'eng': build_qlabel(1),
             'atcrd': build_qlabel('atcrd', 'last'),
             'atnum': build_qlabel('atnum'),
            }
    dfile = DataFile(fname)
    res = {}
    res['fname'] = fname
    try:
        data = dfile.get_data(*dkeys.values())
        atmnum = len(data[dkeys['atnum']]['data'])
        res['apt'] = np.array(data[dkeys['apt']]['data']).reshape(atmnum,3,3)
        res['aat'] = np.array(data[dkeys['aat']]['data']).reshape(atmnum,3,3)
    except:
        data = dfile.get_data(*dkeys2.values())
        
    res['atnum'] = data[dkeys['atnum']]['data']
    res['atlab'] = convert_labsymb(True, *data[dkeys['atnum']]['data'])
    res['atcrd'] = np.array(data[dkeys['atcrd']]['data'])*PHYSFACT.bohr2ang
    res['eng'] = data[dkeys['eng']]['data']
    #res['apt'] = np.array(data[dkeys['apt']]['data']).reshape(len(res['atnum']),3,3)
    #res['aat'] = np.array(data[dkeys['aat']]['data']).reshape(len(res['atnum']),3,3)
    return res

def get_bondsdatatoobg(prefix, suffix, hxobj, nterms):
    res = LocalModes(hxobj.atnum, hxobj.refcrd, nterms)
    bonds = hxobj.hxbonds
    fname = prefix +"_bond_H{}_{:02d}_{:02d}_step"
    fname2 = prefix +"_bond_H{}_{:02d}_{:02d}_step{:03d}_"+suffix+".fchk"
    for bnd in bonds:
        atype = hxobj.getsecatom(bnd)
        lfiles = glob.glob(fname.format(atype, bnd[0]+1, bnd[1]+1)+"*")
        # print(lfiles)
        tmpres = {'eng': [], 'len':[], 'apt1': [], 'aat1': [], 'apt2': [], 'aat2': []}
        for i in range(len(lfiles)):               
            tmp_data =  get_fchk(fname2.format(atype, bnd[0]+1, bnd[1]+1, i))
            #print(tmp_data)
            tmpres['eng'].append(tmp_data['eng'])
            tmpres['len'].append(np.abs(tmp_data['atcrd'][bnd[0], 2]-tmp_data['atcrd'][bnd[1], 2]))
            tmpres['apt1'].append(tmp_data['apt'][bnd[0], 2, :])
            tmpres['aat1'].append(tmp_data['aat'][bnd[0], 2, :])
            tmpres['apt2'].append(tmp_data['apt'][bnd[1], 2, :])
            tmpres['aat2'].append(tmp_data['aat'][bnd[1], 2, :])

        tmp_bond = LmodDeriv()
        tmp_bond.bond = bnd
        tmp_bond.eng = tmpres['eng']
        tmp_bond.blen = tmpres['len']
        tmp_bond.apt = [tmpres['apt1'], tmpres['apt2']]
        tmp_bond.aat = [tmpres['aat1'], tmpres['aat2']]
        res.addbond(tmp_bond)
    return res

