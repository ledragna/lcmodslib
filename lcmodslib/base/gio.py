from estampes.parser import DataFile, build_qlabel
from estampes.tools.atom import convert_labsymb
from estampes.data.physics import PHYSFACT
from estampes.data.atom import atomic_data

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
    template ="""%mem=20GB
%nprocshared=12
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
            fopen.write(template.format(CHK=fout[:-3]+'chk',
                                        FUN=qmdata['functional'],
                                        BASIS=qmdata['basis'],
                                        MOLCHR=qmdata['molchr'],
                                        MOLSPN=qmdata['molspn'],
                                        ADDROOT=qmdata['addroot'],
                                        COMM=tmplcomm.format(gid=i+1, val="test")))
            fopen.write(geomstr)
            fopen.write('\n')
            fopen.write(f'{qmdata["addline"]}\n\n')
        
