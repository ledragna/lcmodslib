"""
Class to manipulate the geometries
"""
import numpy as np

from estampes.tools.atom import convert_labsymb
from estampes.data.atom import atomic_data
from scipy.spatial import distance_matrix

class LocalMolecule():
    """
    Simple class handling the molecular system
    """
    def __init__(self, atnum, refcrd, energy=None, thrs=0.4):
        """
        atnums: list of atomic number
        refcrd: reference atomic crd in Angstrom as natom x 3 matrix
        """
        #self._atn = [element(x).atomic_number for x in atnums]
        if isinstance(atnum[0], int):
            self._atnum = atnum
            self._atlab = convert_labsymb(True, *self._atnum)
        else:
            self._atlab = []
            for x in atnum:
                if len(x) == 1:
                    el = x.upper()
                elif len(x) == 2:
                    el = x[0].upper() + x[1]
                else:
                    print("read element ",x)
                    raise ValueError("I don't know about this one")
                self._atlab.append(el)
            self._atnum = convert_labsymb(False, *self._atnum)
        self._natom = len(self._atnum)
        # BUG think a mre efficent way
        self._refcrd = refcrd
        self._actcrd = refcrd.copy()
        self._thrs = thrs
        # use the covalent radii to evaluate the connectivity
        atdat = atomic_data(*set(self._atlab))
        self._atmass = np.array([atdat[at]['mass'] for at in self._atlab])
        self._ard = np.array([atdat[at]['rcov'][0]/100 for at in self._atlab])
        self._compute_ardmat()
        self._updatecref()
    
    @property
    def actcrd(self):
        return self._actcrd
    
    @property
    def atmass(self):
        return self._atmass
    
    @property
    def atlab(self):
        return self._atlab

    @property
    def atnum(self):
        return self._atnum
    
    @property
    def refcrd(self):
        return self._refcrd
           
    def _compute_ardmat(self):
        """
        compute the sum of covalence radius matrix
        """
        tmp = self._ard[:, np.newaxis] + self._ard[np.newaxis,:]
        self._ardmat = tmp

    def _compute_conn(self, crd):
        """
        evaluate the connectivity of a given coordinate based on the atomic
        number of the reference structure
        crd: coordinates in Angstrom as natom x 3 matrix
        Return:
            A boolen matrix
        """
        dmat = distance_matrix(crd, crd)
        ardsum = self._ardmat + (self._thrs)
        ardsum[np.diag_indices(self._ard.shape[0])] = 0.
        return (ardsum - dmat) > 0

    def _updatecref(self):
        """
        update the connectivity of the reference structure
        """
        self._cref = self._compute_conn(self._refcrd)

    def set_thrs(self, thrs):
        """
        set a different threshold to evaluate connectivity
        """
        self._thrs = thrs
        self._updatecref()


    def get_conmat(self, crd=None):
        """
        Returns the connectivity matrix of a given coordinates.
        If None is given the returns the connectivity matrix of the reference
        """
        if crd is None:
            return self._cref
        return self._compute_conn(crd)
    
    def _indexinthesystem(self, indx):
        if indx < 0 or indx >= self._natom:
            raise ValueError("index out of range")              
    
    def orientzaxisbond(self, bond, atm3=None):
        """
        Orients the z axis along a bond positive towards atm1.
        bond: tuple(atm1,atm2)
        The origin is set in the centre of mass.
        A third atom can be given to define the plane, if None
        get from connectivity
        """
        normy = lambda x: x/np.sqrt(np.dot(x, x))
        # Checks
        self._indexinthesystem(bond[0])
        self._indexinthesystem(bond[1])
        # Checks if bond[1] as connection if not swap bond
        if not (self._cref[bond[1], :] > 0).sum() > 1:
            bond = (bond[1], bond[0])
        # Defines the third atom if not given
        if atm3 is None:
            _tmp = np.where(self._cref[bond[1], :] > 0)[0]
            if _tmp[0] != bond[0]:
                atm3 = _tmp[0]
            else:
                atm3 = _tmp[1]
        else:
            self._indexinthesystem(atm3)
        # compute the centre of mass
        cmass = np.average(self._actcrd[bond, :], axis=0,
                           weights=self._atmass[list(bond)])
        # translat
        _tmpcrd = self._actcrd - cmass
        zaxis = normy(_tmpcrd[bond[0]] - _tmpcrd[bond[1]])
        tplan = normy(_tmpcrd[bond[1]] - _tmpcrd[atm3])
        ytmp = normy(np.cross(zaxis, tplan))
        xtmp = normy(np.cross(ytmp, zaxis))
        rotmat = np.zeros((3,3))
        rotmat[:, 0] = xtmp
        rotmat[:, 1] = ytmp
        rotmat[:, 2] = zaxis
        return _tmpcrd@rotmat
    
    def bondstretch(self, bond, step):
        # check if the molecule is already oriented
        # if (self._actcrd[bond,:2] > 1e-7).any():
        self._actcrd = self.orientzaxisbond(bond)
        mwstep = step * (self._atmass[bond[0]]*self._atmass[bond[1]]/
                         (self._atmass[bond[0]]+self._atmass[bond[1]]))
        # print((1/self._atmass[list(bond)]*mwstep))
        _tmpcrd = self._actcrd.copy()
        _tmpcrd[bond,2] += (1/self._atmass[list(bond)]*np.array([1,-1])*mwstep)
        return _tmpcrd
    
    def redmass(self, bond):
        for i in bond:
            self._indexinthesystem(i)
        return self._atmass[bond[0]]*self._atmass[bond[1]]/(self._atmass[bond[0]]+self._atmass[bond[1]])
        
class XHstreching(LocalMolecule):
    """
    Class to automatically handle the XH stretching within the molecule
    """
    
    def __init__(self, atnum, refcrd, xtype="C", **kwargs):
        super().__init__(atnum, refcrd, **kwargs)
        # checks atmtype
        if not isinstance(xtype, list):
            xtype = [xtype]
        _xtype = []
        for _x in xtype:
            if len(_x) == 1:
                _x = _x.upper()
            elif len(_x) == 2:
                _x = _x[0].upper() + _x[1]
            else:
                print(f"read element {_x}")
                raise ValueError("I don't know about this one")
            if not _x in self._atlab:
                print(f"Element {_x} not in molecular system")
                raise ValueError("Not in molecular system")
            _xtype.append(_x)
        self._xtypes = _xtype
        self._get_ch_bonds()
        
    @property
    def hxbonds(self):
        return self._hxbonds
    
    @property
    def xtypes(self):
        return self._xtypes
    
    def _get_ch_bonds(self):
        xhbonds = []
        xtyplst = []
        hpos = np.where(np.array(self._atnum) == 1)[0]
        for i in hpos:
            _xindx = np.where(self._cref[i, :])[0][0]
            if self._atlab[_xindx] in self._xtypes:
                xhbonds.append((i, _xindx))
                xtyplst.append(self._atlab[_xindx])
        if not xhbonds:
            print(f"No bonds found with the selected atoms")
            raise ValueError("No Bonds selected")
        self._hxbonds = xhbonds
        self._xtyplst = xtyplst
    
    def hxstretch(self, bndindx, step):
        """
        stretch a HX based on internal bndx indexing
        """
        _bond = self._checkbond(bndindx)
        return self.bondstretch(self._hxbonds[_bond], step)
    
    def _checkbond(self, bond):
        if isinstance(bond, tuple):
            if bond in self._hxbonds:
                try:
                    _bond = self._hxbonds.index(bond)
                except ValueError:
                    try:
                        _bond = self._hxbonds.index((bond[1],bond[0]))
                    except ValueError:
                        raise ValueError(f"{bond} not in the list of selected bonds")
        else:
            _bond = int(bond)
            if not (_bond > 0 and _bond < len(self._hxbonds)):
                raise ValueError("index out of range")
        return _bond              
    
    def scanhx(self, bond, lower=-0.33, step=0.017, nstep=50):
        """
        bond the bond to stretch can be the inter index or a tuple
        lower: lower bound of the scan
        step: the step size
        nstep: the number of step
        """
        # in order to get the equilibrium position the lower bounds is corrected
        _lower = (np.abs(lower)//step)*step*lower/np.abs(lower)
        # check bond
        _bond = self._checkbond(bond)
        geoms = []
        for i in range(nstep):
            xstep = i*step + _lower
            geoms.append(self.bondstretch(self._hxbonds[_bond], xstep))
        
        return geoms

    def getsecatom(self, bond=None):
        if bond is None:
            return self._xtyplst
        _bond = self._checkbond(bond)
        return self._xtyplst[_bond]
        
