import numpy as np
import typing as tp

from estampes.tools.atom import convert_labsymb
from estampes.data.atom import atomic_data
from scipy.spatial import distance_matrix

def readlistnum(x: str, nmax: tp.Optional[tp.Union[None, int]]=None) -> tp.List[int]:
    result = []
    for part in x.split(','):
        if '-' in part:
            a, b = part.split('-')
            if a == "":
                a = 0
            else:
                a = int(a)
            if b == "":
                if nmax:
                    b = int(nmax)
                else:
                    raise NotImplementedError()
            else:
                b = int(b)
            result.extend(range(a, b + 1))
        else:
            a = int(part)
            result.append(a)
    return result

def connectivitymatrix(atnum: tp.Union[np.ndarray, list], refcrd: np.ndarray, thrs: float=0.4) -> np.ndarray:
    if isinstance(atnum[0], int):
            _atnum = atnum
            _atlab = convert_labsymb(True, *_atnum)
    else:
        _atlab = []
        for x in atnum:
            if len(x) == 1:
                el = x.upper()
            elif len(x) == 2:
                el = x[0].upper() + x[1]
            else:
                print("read element ",x)
                raise ValueError("I don't know about this one")
            _atlab.append(el)
        _atnum = convert_labsymb(False, * atnum)
    _natom = len(_atnum)
    # BUG think a mre efficent way
    _refcrd = refcrd
    _thrs = thrs
    # use the covalent radii to evaluate the connectivity
    atdat = atomic_data(*set(_atlab))
    # _atmass = np.array([atdat[at]['mass'] for at in _atlab])
    _ard = np.array([atdat[at]['rcov'][0]/100 for at in _atlab])
    _ardmat = _ard[:, np.newaxis] + _ard[np.newaxis,:]
    dmat = distance_matrix(_refcrd, _refcrd)
    ardsum = _ardmat + _thrs
    ardsum[np.diag_indices(_natom)] = 0.
    return (ardsum - dmat) > 0

def angle(rvec1: np.ndarray, rvec2: np.ndarray) -> float:
    cos_alpha = np.dot(rvec1, rvec2)
    sin_alpha = np.linalg.norm(np.cross(rvec1, rvec2))
    alpha = np.arctan2(sin_alpha, cos_alpha)
    #return np.degrees(alpha)
    return np.rad2deg(alpha)
