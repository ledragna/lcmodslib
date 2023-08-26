import numpy as np
from estampes.data.atom import atomic_data
from estampes.tools.atom import convert_labsymb
from lcmodslib.base.utils import connectivitymatrix, angle

class LmodDeriv():
    """A class to store the data from QM calculations of a defined modes
    """
    
    def __init__(self):
        self._eng = None
        self._blen = None
        self._bond = None
        self._apt = None
        self._aptmask = None
        self._aat = None
        self._aatmask = None
        self._ams = None
        self._rmas = None
        self._rdis = None
        
    @property
    def eng(self):
        return self._eng
    
    @property
    def blen(self):
        return self._blen
    
    @property
    def apt(self):
        return self._apt
    
    @property
    def aat(self):
        return self._aat
    
    @property
    def bond(self):
        return self._bond
    
    @property
    def ams(self):
        return self._ams
    
    @property
    def rdis(self):
        return self._rdis
            
    @eng.setter
    def eng(self, val):
        self._eng = val

    @blen.setter
    def blen(self, val):
        self._blen = val
        
    @apt.setter
    def apt(self, val):
        _tmp_mask, _tmp_apt = self._checktensor(val)
        self._apt = _tmp_apt
        self._aptmask = _tmp_mask

    @aat.setter
    def aat(self, val):
        _tmp_mask, _tmp_aat = self._checktensor(val)
        self._aat = _tmp_aat
        self._aatmask = _tmp_mask
    
    @bond.setter
    def bond(self, val):
        self._bond = val
    
    @ams.setter
    def ams(self, val):
        self._rmas = (val[0]*val[1])/(val[0]+val[1])
        self._ams = val
    
    @rdis.setter
    def rdis(self, val):
        self._rdis = val

    @staticmethod
    def _checktensor(val):
        _mask = []
        _tmp = [[], []]
        for k, nvl in enumerate(val[0]):
            if not nvl is None:
                _mask.append(k)
                _tmp[0].append(nvl)
                _tmp[1].append(val[1][k])
        return (_mask, [np.array(_tmp[0]), np.array(_tmp[1])]) 

    def lmanharm(self, deg: int=8):
        slight = 2.99792458e10 # cm/s
        avo = 6.02214076e23 # mol -1
        plank = 6.62607015e-27 # erg*s
        conv = 1e16
        hartree2erg = 4.3597447222071e-11
        # checks if equilibrium is in the scan
        if not self._rdis is None:
            _rvec = self.blen - self._rdis
        else:
            minpos = self.blen[np.array(self.eng).argmin()]  
            _rvec = self.blen - minpos
        poli = np.polynomial.polynomial.polyfit(_rvec, self.eng, deg)
        k2 = poli[2]*2
        k3 = poli[3]*6
        k4 = poli[4]*24
        #print(poli)
        #print(k2,k3,k4)
        # X = h/(64pi^2*mr*c)*(5/3Kiii^2/Kii^2-Kiiii/Kii)
        # X = hbar/(32pi*mr*c)*(5/3Kiii^2/Kii^2-Kiiii/Kii)
        # w = 1/(2pi*c)*sqrt(Kii/mr)
        wau = 1/(2*np.pi*slight)*np.sqrt(k2*hartree2erg*conv/(self._rmas/avo))
        xau = plank*conv/(64*np.pi**2*self._rmas/avo*slight)*(5/3*k3**2/k2**2-k4/k2)
        return (wau, xau)
    
    @staticmethod
    def compute_integrals(omega: float, chi: float, rdmass: float, quanta: int) -> list[list[float]]:
        """Computes the transition <0|z|n>,<0|z^2|n>,<0|z^3|n> and <0|p|0>,<0|zp|n>,<0|z^2p|n> transition integrals see:

        Args:
            omega (float): w_0 from morse potential fitting (cm-1)
            chi (float): X from morse potential fitting (cm-1)
            rdmass (float): bi-atomic systems reduced mass (amu)
            quanta (int): quanta of interest (1 to 6)

        Returns:
            list[list[float]]: [[q,q^2,q^3],[p, qp, q^2p]]
        """
        integrals = {1:{}, 2:{}, 3:{}, 4:{}, 5:{}, 6:{}}
        # Adding q4 and q3p harmonic for test
        # 1
        integrals[1]['q'] = lambda d, k: d/(2*np.pi)*(1+k/2)
        integrals[1]['q2'] = lambda d, k: (5/4)*d**2/np.pi**2*np.sqrt(k)*(1+52/30*k)
        integrals[1]['q3'] = lambda d, k: (3/8)*d**3/np.pi**3*(1+37/4*k)
        # integrals[1]['q4'] = lambda d, k: 0
        integrals[1]['p'] = lambda d, k: -hbar*np.pi/d*(1-3/2*k)
        integrals[1]['qp'] = lambda d, k: -hbar*5/4*np.sqrt(k)*(1-4/15*k)
        integrals[1]['q2p'] = lambda d, k: hbar*d/(4*np.pi)*(1-25/4*k)
        # integrals[1]['q3p'] = lambda d, k: 0
        # 2
        integrals[2]['q'] = lambda d, k: -d/(2*np.sqrt(2)*np.pi)*np.sqrt(k)*(1+3/2*k)
        integrals[2]['q2'] = lambda d, k: d**2/(2*np.sqrt(2)*np.pi**2)*(1-2*k)
        integrals[2]['q3'] = lambda d, k: 9/4*d**3/(np.sqrt(2)*np.pi**2)*np.sqrt(k)*(1+65/72*k)
        # integrals[2]['q4'] = lambda d, k: (3*d**4)/(2**(5/2)*np.pi**4) 
        integrals[2]['p'] = lambda d, k: hbar*np.sqrt(2)*np.pi/d*np.sqrt(k)*(1-3/2*k)
        integrals[2]['qp'] = lambda d, k: -hbar*np.sqrt(2)/2*(1-5*k)
        integrals[2]['q2p'] = lambda d, k: -hbar*7*d/(2*np.sqrt(2)*np.pi)*np.sqrt(k)*(1-19/12*k)
        # integrals[2]['q3p'] = lambda d, k: 0
        # 3
        integrals[3]['q'] = lambda d, k: np.sqrt(2/3)*d/(2*np.pi)*k*(1+3*k)
        integrals[3]['q2'] = lambda d, k: -np.sqrt(3/2)*d**2/(2*np.pi**2)*np.sqrt(k)
        integrals[3]['q3'] = lambda d, k: np.sqrt(3/2)*d**3/(4*np.pi**3)*(1-19/2*k)
        # integrals[3]['q4'] = lambda d, k: 0
        integrals[3]['p'] = lambda d, k: -hbar*np.sqrt(3/2)*2*np.pi/d*k*(1-k)
        integrals[3]['qp'] = lambda d, k: hbar*np.sqrt(3/2)*3/2*np.sqrt(k)*(1-4*k)
        integrals[3]['q2p'] = lambda d, k: -hbar*np.sqrt(3/2)*d/(2*np.pi)*(1-85/6*k)
        # integrals[3]['q3p'] = lambda d, k: 0
        # 4
        integrals[4]['q'] = lambda d, k: -np.sqrt(3/2)*d/(2*np.pi)*k**(3/2)*(1+5*k)
        integrals[4]['q2'] = lambda d, k: 11*d**2/(4*np.sqrt(6)*np.pi**2)*k*(1+2*k)
        integrals[4]['q3'] = lambda d, k: -3*np.sqrt(3/2)*d**3/(4*np.pi**3)*np.sqrt(k)*(1-55/12*k)
        # integrals[4]['q4'] = lambda d, k: (np.sqrt(6)*d**4)/(8*np.pi**4)
        integrals[4]['p'] = lambda d, k: hbar*2*np.sqrt(6)*np.pi/d*k**(3/2)
        integrals[4]['qp'] = lambda d, k: -hbar*11/np.sqrt(6)*k*(1-3*k)
        integrals[4]['q2p'] = lambda d, k: hbar * np.sqrt(6)*d/np.pi *np.sqrt(k) *(1-59/6*k)
        # integrals[4]['q3p'] = lambda d, k: 0
        # 5
        integrals[5]['q'] = lambda d, k: np.sqrt(6/5)*d/np.pi*k**2*(1+15/2*k)
        integrals[5]['q2'] = lambda d, k: -5*np.sqrt(5/6)*d**2/(2*np.pi**2)*k**(3/2)*(1+219/50*k)
        integrals[5]['q3'] = lambda d, k: 7*np.sqrt(15/2)*d**3/(8*np.pi**3)*k*(1-101/70*k)
        # integrals[5]['q4'] = lambda d, k: 0
        integrals[5]['p'] = lambda d, k: -hbar*2*np.sqrt(30)*np.pi/d*k**2
        integrals[5]['qp'] = lambda d, k: hbar*25/2*np.sqrt(5/6)*k**(3/2)*(1-81/50*k)
        integrals[5]['q2p'] = lambda d, k: -hbar*35/4*np.sqrt(5/6)*d/np.pi*k*(1-379/50*k)
        # integrals[5]['q3p'] = lambda d, k: 0
        # 6
        integrals[6]['q'] = lambda d, k: -np.sqrt(5)*d/np.pi*k**(5/2)
        integrals[6]['q2'] = lambda d, k: 137*d**2/(12*np.sqrt(5)*np.pi**2)*k**2
        integrals[6]['q3'] = lambda d, k: -45*np.sqrt(5)*d**3/(16*np.pi**3)*k**(3/2)*(1+147/90*k)
        # integrals[6]['q4'] = lambda d, k: 0
        integrals[6]['p'] = lambda d, k: hbar*12*np.sqrt(5)*np.pi/d*k**(5/2)
        integrals[6]['qp'] = lambda d, k: -hbar*137*np.sqrt(5)/10*k**2
        integrals[6]['q2p'] = lambda d, k: hbar*45*np.sqrt(5)*d/(np.pi*4)*k**(3/2)
        # integrals[6]['q3p'] = lambda d, k: 0

        hbar = 1.054571817e-27 # erg*s
        planck = 6.62607015e-27 # erg*s
        slight = 2.99792458e10 # cm/s
        # avo = 6.02214076e23 # mol -1
        da2gra = 1.66053906660e-24 # g
        cm2ang = 1e8
        dq = planck/(2*slight*rdmass*da2gra)/omega
        _k = chi/omega
        # d in angstrom
        _d = np.sqrt(dq)*cm2ang
    
        return [[integrals[quanta]['q'](_d, _k), integrals[quanta]['q2'](_d, _k), integrals[quanta]['q3'](_d, _k)],
                [-integrals[quanta]['p'](_d, _k), -integrals[quanta]['qp'](_d, _k), -integrals[quanta]['q2p'](_d, _k)]]    

    def _compute_dipoles(self, quanta, omega=None, chi=None,
                         deg=8, terms=3):
        """
        First step interpolation for each atom and each component
        si tiene 0 1 2
        
        d^2 = h/(2*c*mr)* 1/ωl
        
        """
        e2esu = 4.80320427e-10    
        # AAT
        # bohr2cm = 5.2917721090e-9
        bohrmagneton = 9.274010078e-21 # egr/G
        hbar = 1.054571817e-27
        aat_c = 2*bohrmagneton/e2esu**2*hbar
        ang2cm = 1e-8
        da2gra = 1.66053906660e-24
        _inter_apt = np.zeros((2, 3, 3))
        _inter_aat = np.zeros((2, 3, 3))
        if not self._rdis is None:
            _rvec = self.blen - self._rdis
        else:
            # Possible BUG check me
            minpos = self.beln[np.array(self.eng).argmin()]  
            _rvec = self.blen - minpos
        ## Expansion coefficients to be check!!!!
        expcoef = np.array([[1,1/2,1/6],
        #expcoef = np.array([[1,1,1/2],
                            [1,1,1/2]])
        # print(np.array(apt[0]))
        for i in range(2):
            for j in range(3):
                _tmp = np.polynomial.polynomial.polyfit(_rvec[self._aptmask], self.apt[i][:, j], deg)
                _inter_apt[i,j,:] = _tmp[:3]
                _tmp = np.polynomial.polynomial.polyfit(_rvec[self._aatmask], self.aat[i][:, j], deg)
                _inter_aat[i,j,:] = _tmp[:3]
    
        tz = np.array([self._ams[1]/(self._ams[0]+self._ams[1]), -self._ams[0]/(self._ams[0]+self._ams[1])])
        # tz / mr 
        omass = np.array([1/self._ams[0], -1/self._ams[1]])*1/da2gra
        if omega == None or chi == None:
            tmp = self.lmanharm(deg=deg)
            omega = tmp[0]
            chi = tmp[1]
        intgs = np.array(self.compute_integrals(omega, chi, self._rmas, quanta))
        intgs *= expcoef
        if terms - 1 < 2:
            intgs[:, terms: ] = 0
        # APTS
        _inter_apt *= tz[:, np.newaxis, np.newaxis]
        _inter_apt = _inter_apt.sum(axis=0)        
        _inter_apt *= intgs[0,:][np.newaxis, :]
        #_mu = _inter_apt.sum(axis=1)
        # AAT
        _inter_aat *= omass[:, np.newaxis, np.newaxis]
        _inter_aat = _inter_aat.sum(axis=0)
        _inter_aat *= intgs[1,:][np.newaxis, :]
        #_m = _inter_aat.sum(axis=1)
        return(_inter_apt*e2esu*ang2cm, _inter_aat*aat_c/ang2cm)
        # return (_mu*e2esu*ang2cm, _m*aat_c/ang2cm)

    def compute_intensities(self, quanta, omega=None, chi=None,
                            deg=8, terms=3, debug=False):
        """
        First step interpolation for each atom and each component
        si tiene 0 1 2
        
        d^2 = h/(2*c*mr)* 1/ωl
        
        """
        _inter_apt, _inter_aat = self._compute_dipoles(quanta, omega, chi, deg,
                                                       terms)
        # APT
        _mu = _inter_apt.sum(axis=1)
        # AAT
        _m = _inter_aat.sum(axis=1)
        return (_mu, _m)

    
class LocalModes():
    """Class to handle the local modes of a molecular system
    """
    def __init__(self, atnum, refcrd, terms=3):
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
        self._cref = connectivitymatrix(self._atnum,
                                        self._refcrd)
        # use the covalent radii to evaluate the connectivity
        atdat = atomic_data(*set(self._atlab))
        self._atmass = np.array([atdat[at]['mass'] for at in self._atlab])
        # Number of terms
        self._ntr = terms
        self._bnd = []
        self._drv = {}
        self._res = {}
    
    def _indexinthesystem(self, indx):
        if indx < 0 or indx >= self._natom:
            raise ValueError("index out of range")
            
    def _refdistance(self, bond):
        rvec = self._refcrd[bond[0]] - self._refcrd[bond[1]]
        return np.sqrt(np.dot(rvec, rvec))
            
    def addbond(self, bond, deg=8):
        bndindx = bond.bond
        atmass = []
        for i in bndindx:
            self._indexinthesystem(i)
            atmass.append(self._atmass[i])
        refdist = self._refdistance(bndindx)
        if bond.ams is None:
            bond.ams = atmass
        if bond.rdis is None:
            minpos = bond.blen[np.array(bond.eng).argmin()] 
            if np.abs(minpos - refdist) < 1e-6:
                bond.rdis = minpos
            else: 
                bond.rdis = refdist
        self._bnd.append(bndindx)
        self._drv[bndindx] = bond
        self._res[bndindx] = {}
        omega, chi = bond.lmanharm()
        # print(omega, chi)
        self._res[bndindx]['omega'] = omega # [0]
        self._res[bndindx]['chi'] = chi # [0]
        # self._res[bndindx]['frq'] = [omega[0]-2*chi[0], 2*omega[0]-6*chi[0], 3*omega[0]-12*chi[0]]
        self._res[bndindx]['frq'] = [i*omega-(i**2+i)*chi for i in range(1,7)]
        self._res[bndindx]['ds'] = []
        self._res[bndindx]['rs'] = []
        self._res[bndindx]['aa'] = []
        self._res[bndindx]['m2'] = []
        self._res[bndindx]['hds'] = []
        self._res[bndindx]['hrs'] = []
        self._res[bndindx]['haa'] = []
        self._res[bndindx]['hm2'] = []
        # BUG 
        try:
            for i in range(6):
                tmpdip = bond.compute_intensities(i+1, omega=omega, chi=chi, terms=self._ntr)
                # Harmonic integrals only
                tmphdip = bond.compute_intensities(i+1, omega=omega, chi=0, terms=self._ntr)
                self._res[bndindx]['ds'].append(np.dot(tmpdip[0], tmpdip[0]))
                self._res[bndindx]['rs'].append(np.dot(tmpdip[0], tmpdip[1]))
                self._res[bndindx]['aa'].append(angle(tmpdip[0], tmpdip[1]))
                self._res[bndindx]['m2'].append(np.dot(tmpdip[1], tmpdip[1]))
                self._res[bndindx]['hds'].append(np.dot(tmphdip[0], tmphdip[0]))
                self._res[bndindx]['hrs'].append(np.dot(tmphdip[0], tmphdip[1]))
                self._res[bndindx]['haa'].append(angle(tmphdip[0], tmphdip[1]))
                self._res[bndindx]['hm2'].append(np.dot(tmphdip[1], tmphdip[1]))
        except IndexError:
            pass
            
    def printlmodes(self, quanta):
        """
        """
        print(16*"#")
        print("{:^16s}".format(f"nu = {quanta}"))
        print(16*"#")
        print("{:^8s}{:^9s}{:^12s}{:^12s}{:^9s}".format("Bond", "Freq.", "DS", "RS", "angle"))
        line ="{a[0]:^4d}{a[1]:^4d}{b:9.2f}{c:12.4E}{d:12.4E}{e:9.2f}"
        for bnx in self._bnd:
            print(line.format(a=bnx, b=self._res[bnx]['frq'][quanta-1],
                              c=self._res[bnx]['ds'][quanta-1],
                              d=self._res[bnx]['rs'][quanta-1],
                              e=self._res[bnx]['aa'][quanta-1]))


    def lmodes2string(self, maxquanta, head=True, harmdip=False):
        """
        """
        keys = ['ds', 'rs', 'aa', 'm2']
        if harmdip:
            keys = ['h'+x for x in keys]
        string = ""
        if head:
            string += "#{:^6s}{:^14s}{:^9s}{:^12s}{:^12s}{:^9s}{:^12s}\n".format("Qnt","Bond", "Freq.", "DS", "RS", "angle", "m dot m")
        line =" {e:^6d}{a[0]:^4d}{f[0]:^3s}{a[1]:^4d}{f[1]:^3s}{b:9.2f}{c:12.4E}{d:12.4E}{g:9.2f}{h:12.4E}\n"
        for qnt in range(maxquanta):
            for bnx in self._bnd:
                bnxt = [x+1 for x in list(bnx)]
                _atlb = [self._atlab[x] for x in list(bnx)]
                string += line.format(a=bnxt,
                                      b=self._res[bnx]['frq'][qnt],
                                      c=self._res[bnx][keys[0]][qnt],
                                      d=self._res[bnx][keys[1]][qnt],
                                      e=qnt+1, f=_atlb,
                                      g=self._res[bnx][keys[2]][qnt],
                                      h=self._res[bnx][keys[3]][qnt])
        return string

    def omgchi2string(self, head=True):
        """
        
        """
        string = ""
        if head:
            string += "#{:^14s}{:^9s}{:^9s}\n".format("Bond", "Omega", "Chi")
        line = " {a[0]:^4d}{d[0]:^3s}{a[1]:^4d}{d[1]:^3s}{b:9.2f}{c:9.2f}\n"
        for bnx in self._bnd:
            bnxt = [x+1 for x in list(bnx)]
            _atlb = [self._atlab[x] for x in list(bnx)]
            string += line.format(a=bnxt,
                                  b=self._res[bnx]['omega'],
                                  c=self._res[bnx]['chi'],
                                  d=_atlb)
        return string

 
    
    def get_lmodes(self, spctyp, quanta):
        eng = []
        spc = []
        for bnx in self._bnd:
            eng.append(self._res[bnx]['frq'][quanta-1])
            spc.append(self._res[bnx][spctyp][quanta-1])
        return [eng, spc]
    
    def getlocal2wrld(self, bond, atm3=None):
        """_summary_

        Args:
            bond (_type_): _description_
            atm3 (_type_, optional): _description_. Defaults to None.
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
        wrld2local = np.identity(4)
        wrld2local[:3, 0] = xtmp
        wrld2local[:3, 1] = ytmp
        wrld2local[:3, 2] = zaxis
        wrld2local[:3, 3] = cmass
        # print(wrld2local)
        return np.linalg.inv(wrld2local)
 

class LocalConfs():
    """
    Class to handle multiple conformers
    """

    def __init__(self, natms, bonds) -> None:
        self._natms = natms
        self._bnds = bonds
        self._nconf = 0
        self._engs = []

    def addconf(self, conf):
        if conf._natom != self._natms:
            print("Different molecules")
        elif self._bnds != conf._bnd:
            print("Differnt bonds")
        else:
            pass


    