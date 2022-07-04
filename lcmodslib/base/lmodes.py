import numpy as np

class LmodDeriv():
    
    def __init__(self):
        self._eng = None
        self._blen = None
        self._bond = None
        self._apt = None
        self._aat = None
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
        self._apt = val

    @aat.setter
    def aat(self, val):
        self._aat = val
    
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
    
    def lmanharm(self, deg=8):
        slight = 2.99792458e10 # cm/s
        avo = 6.02214076e23 # mol -1
        plank = 6.62607015e-27 # erg*s
        conv = 1e16
        hartree2erg = 4.3597447222071e-11
        # checks if equilibrium is in the scan
        if not self._rdis is None:
            _rvec = self.blen - self._rdis
        else:
            minpos = self.beln[np.array(self.eng).argmin()]  
            _rvec = selb.blen - minpos
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
    def compute_integrals(omega, chi, rdmass, quanta):
        integrals = {1:{}, 2:{}, 3:{}}
        # 1
        integrals[1]['q'] = lambda d, k: d/(2*np.pi)*(1+k/2)
        integrals[1]['q2'] = lambda d, k: (5/4)*d**2/np.pi**2*np.sqrt(k)*(1+52/30*k)
        integrals[1]['q3'] = lambda d, k: (3/8)*d**3/np.pi**3*(1+37/4*k)
        integrals[1]['p'] = lambda d, k: -hbar*np.pi/d*(1-3/2*k)
        integrals[1]['qp'] = lambda d, k: -hbar*5/4*np.sqrt(k)*(1-4/15*k)
        # only harmonic
        integrals[1]['q2p'] = lambda d, k: hbar*d/4
        # 2
        integrals[2]['q'] = lambda d, k: -d/(2*np.sqrt(2)*np.pi)*np.sqrt(k)*(1+3/2*k)
        integrals[2]['q2'] = lambda d, k: d**2/(2*np.sqrt(2)*np.pi**2)*(1-2*k)
        integrals[2]['q3'] = lambda d, k: 9/4*d**3/(np.sqrt(2)*np.pi**2)*np.sqrt(k)*(1+65/72*k)
        integrals[2]['p'] = lambda d, k: hbar*np.sqrt(2)*np.pi/d*np.sqrt(k)*(1-3/2*k)
        integrals[2]['qp'] = lambda d, k: -hbar*np.sqrt(2)/2*(1-5*k)
        # only harmonic
        integrals[2]['q2p'] = lambda d, k: 0
        # 3
        integrals[3]['q'] = lambda d, k: np.sqrt(2/3)*d/(2*np.pi)*k*(1+3*k)
        integrals[3]['q2'] = lambda d, k: -np.sqrt(3/2)*d**2/(2*np.pi**2)*np.sqrt(k)
        integrals[3]['q3'] = lambda d, k: np.sqrt(3/2)*d**3/(4*np.pi**3)*(1+31/2*k)
        integrals[3]['p'] = lambda d, k: -hbar*np.sqrt(3/2)*2*np.pi/d*k*(1-k)
        integrals[3]['qp'] = lambda d, k: hbar*np.sqrt(3/2)*3/2*np.sqrt(k)*(1-4*k)
        # only harmonic
        integrals[3]['q2p'] = lambda d, k: hbar*np.sqrt(3/2)*d/2
        
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

    def conpute_intensities(self, quanta, omega=None, chi=None, deg=8, terms=3):
        """
        First step interpolation for each atom and each component
        si tiene 0 1 2
        
        d^2 = h/(2*c*mr)* 1/ωl
        
        """
        e2esu = 4.80320427e-10    
        # AAT
        # 
        bohr2cm = 5.2917721090e-9
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
            minpos = self.beln[np.array(self.eng).argmin()]  
            _rvec = selb.blen - minpos
        ## Expansion coefficents to be check!!!!
        expcoef = np.array([[1,1/2,1/6],
        #expcoef = np.array([[1,1,1/2],
                            [1,1,1/2]])
        # print(np.array(apt[0]))
        for i in range(2):
            for j in range(3):
                _tmp = np.polynomial.polynomial.polyfit(_rvec, np.array(self.apt[i])[:, j], deg)
                _inter_apt[i,j,:] = _tmp[:3]
                _tmp = np.polynomial.polynomial.polyfit(_rvec, np.array(self.aat[i])[:, j], deg)
                _inter_aat[i,j,:] = _tmp[:3]
    
        tz = np.array([self._ams[1]/(self._ams[0]+self._ams[1]), -self._ams[0]/(self._ams[0]+self._ams[1])])
        omass = np.array([1/self._ams[0], -1/self._ams[1]])*1/da2gra
        if omega == None or chi == None:
            tmp = self.lmanharm(deg=deg)
            omega = tmp[0][0]
            chi = tmp[1][0]
        intgs = np.array(self.compute_integrals(omega, chi, self._rmas, quanta))
        intgs *= expcoef
        if terms - 1 < 2:
            intgs[:, terms: ] = 0
        # APTS
        _inter_apt *= tz[:, np.newaxis, np.newaxis]
        _inter_apt = _inter_apt.sum(axis=0)
        
        _inter_apt *= intgs[0,:][np.newaxis, :]
        # print(_inter_apt)
        _mu = _inter_apt.sum(axis=1)
        # AAT
        _inter_aat *= omass[:, np.newaxis, np.newaxis]
        _inter_aat = _inter_aat.sum(axis=0)
        _inter_aat *= intgs[1,:][np.newaxis, :]
        _m = _inter_aat.sum(axis=1)
        return (_mu*e2esu*ang2cm, _m*aat_c/ang2cm)

    
class LocalModes():
    def __init__(self, atnum, refcrd, terms=3):
        if isinstance(atnum[0], int):
            self._atnum = atnum
            self._atlab = convert_labsymb(True, *self._atnum)
        else:
            self._atlab = []
            for x in atnums:
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
        self._res[bndindx]['omega'] = omega[0]
        self._res[bndindx]['chi'] = chi[0]
        self._res[bndindx]['frq'] = [omega[0]-2*chi[0], 2*omega[0]-6*chi[0], 3*omega[0]-12*chi[0]]
        self._res[bndindx]['ds'] = []
        self._res[bndindx]['rs'] = []
        for i in range(3):
            tmpdip = bond.conpute_intensities(i+1, omega=omega[0], chi=chi[0], terms=self._ntr)
            self._res[bndindx]['ds'].append(np.dot(tmpdip[0], tmpdip[0]))
            self._res[bndindx]['rs'].append(np.dot(tmpdip[0], tmpdip[1]))
            
    def printlmodes(self, quanta):
        """
        """
        print(16*"#")
        print("{:^16s}".format(f"nu = {quanta}"))
        print(16*"#")
        print("{:^8s}{:^9s}{:^12s}{:^12s}".format("Bond", "Freq.", "DS", "RS"))
        line ="{a[0]:^4d}{a[1]:^4d}{b:9.2f}{c:12.4E}{d:12.4E}"
        for bnx in self._bnd:
            print(line.format(a=bnx, b=self._res[bnx]['frq'][quanta-1],
                              c=self._res[bnx]['ds'][quanta-1],
                              d=self._res[bnx]['rs'][quanta-1]))
    
    def get_lmodes(self, spctyp, quanta):
        eng = []
        spc = []
        for bnx in self._bnd:
            eng.append(self._res[bnx]['frq'][quanta-1])
            spc.append(self._res[bnx][spctyp][quanta-1])
        return [eng, spc]